/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/

#include <assert.h>

#if defined(_OPENMP)
#include <omp.h>
#define OMP(x) OMP_(x)
#define OMP_(x) _Pragma(x)
#else
#define OMP(x)
static inline int omp_get_num_threads (void) { return 1; }
static inline int omp_get_thread_num (void) { return 0; }
#endif

#include "graph.h"
#include "../config.h"
#include "../util/arr.h"
#include "../util/qsort.h"
#include "../GraphBLASExt/GxB_Delete.h"
#include "../util/rmalloc.h"
#include "../util/datablock/oo_datablock.h"

#if !defined(NDEBUG)
#define INFO_SUCCESS(info) do { if (info != GrB_SUCCESS) fprintf (stderr, "%s:%d: info: %d error %s\n", __FILE__, __LINE__, (int)info, GrB_error ()); assert(info == GrB_SUCCESS); } while (0)
#else
#define INFO_SUCCESS(info)
#endif

// GraphBLAS Select operator to free edge arrays and delete edges.
static GxB_SelectOp _select_delete_edges = NULL;

/* ========================= Forward declarations  ========================= */
void _MatrixResizeToCapacity(const Graph *g, RG_Matrix m);


/* ========================= RG_Matrix functions =============================== */

// Creates a new matrix;
static RG_Matrix RG_Matrix_New(GrB_Type data_type, GrB_Index nrows, GrB_Index ncols) {
	RG_Matrix matrix = rm_calloc(1, sizeof(_RG_Matrix));
	GrB_Info matrix_res = GrB_Matrix_new(&matrix->grb_matrix, data_type, nrows, ncols);
	assert(matrix_res == GrB_SUCCESS);
	int mutex_res = pthread_mutex_init(&matrix->mutex, NULL);
	assert(mutex_res == 0);
	return matrix;
}

// Returns underlying GraphBLAS matrix.
static inline GrB_Matrix RG_Matrix_Get_GrB_Matrix(RG_Matrix matrix) {
	return matrix->grb_matrix;
}

// Locks the matrix.
static inline void RG_Matrix_Lock(RG_Matrix matrix) {
	pthread_mutex_lock(&matrix->mutex);
}

// Unlocks the matrix.
static inline void _RG_Matrix_Unlock(RG_Matrix matrix) {
	pthread_mutex_unlock(&matrix->mutex);
}

// Free RG_Matrix.
static void RG_Matrix_Free(RG_Matrix matrix) {
	GrB_Matrix_free(&matrix->grb_matrix);
	pthread_mutex_destroy(&matrix->mutex);
	rm_free(matrix);
}

/* ========================= Synchronization functions ========================= */

/* Acquire a lock that does not restrict access from additional reader threads */
void Graph_AcquireReadLock(Graph *g) {
	pthread_rwlock_rdlock(&g->_rwlock);
}

/* Acquire a lock for exclusive access to this graph's data */
void Graph_AcquireWriteLock(Graph *g) {
	pthread_rwlock_wrlock(&g->_rwlock);
	g->_writelocked = true;
}

/* Release the held lock */
void Graph_ReleaseLock(Graph *g) {
	/* Set _writelocked to false BEFORE unlocking
	 * if this is a reader thread no harm done,
	 * if this is a writer thread the writer is about to unlock so once again
	 * no harm done, if we set `_writelocked` to false after unlocking it is possible
	 * for a reader thread to be considered as writer, performing illegal access to
	 * underline matrices, consider a context switch after unlocking `_rwlock` but
	 * before setting `_writelocked` to false. */
	g->_writelocked = false;
	pthread_rwlock_unlock(&g->_rwlock);
}

/* Writer request access to graph. */
void Graph_WriterEnter(Graph *g) {
	pthread_mutex_lock(&g->_writers_mutex);
}

/* Writer release access to graph. */
void Graph_WriterLeave(Graph *g) {
	pthread_mutex_unlock(&g->_writers_mutex);
}

/* Force execution of all pending operations on a matrix. */
static inline void _Graph_ApplyPending(GrB_Matrix m) {
        GrB_Info info;
        info = GrB_wait ();
	INFO_SUCCESS(info);
}

/* ========================= Graph utility functions ========================= */

// Return number of nodes graph can contain.
size_t _Graph_NodeCap(const Graph *g) {
	return g->nodes->itemCap;
}

// Return number of nodes graph can contain.
size_t _Graph_EdgeCap(const Graph *g) {
	return g->edges->itemCap;
}

// Locates edges connecting src to destination.
void _Graph_GetEdgesConnectingNodes(const Graph *g, NodeID src, NodeID dest, int r, Edge **edges) {
	assert(g && src < Graph_RequiredMatrixDim(g) && dest < Graph_RequiredMatrixDim(g) &&
		   r < Graph_RelationTypeCount(g));

	Edge e;
	EdgeID edgeId;
	e.relationID = r;
	e.srcNodeID = src;
	e.destNodeID = dest;

	// relation map, maps (src, dest, r) to edge IDs.
	GrB_Matrix relation = Graph_GetRelationMatrix(g, r);
	GrB_Info res = GrB_Matrix_extractElement_UINT64(&edgeId, relation, src, dest);

	// No entry at [dest, src], src is not connected to dest with relation R.
	if(res == GrB_NO_VALUE) return;

	if(SINGLE_EDGE(edgeId)) {
		// Discard most significate bit.
		edgeId = SINGLE_EDGE_ID(edgeId);
		e.entity = DataBlock_GetItem(g->edges, edgeId);
		e.id = edgeId;
		assert(e.entity);
		*edges = array_append(*edges, e);
	} else {
		/* Multiple edges connecting src to dest,
		 * entry is a pointer to an array of edge IDs. */
		EdgeID *edgeIds = (EdgeID *)edgeId;
		int edgeCount = array_len(edgeIds);

		for(int i = 0; i < edgeCount; i++) {
			edgeId = edgeIds[i];
			e.entity = DataBlock_GetItem(g->edges, edgeId);
			e.id = edgeId;
			assert(e.entity);
			*edges = array_append(*edges, e);
		}
	}
}

// Tests if there's an edge of type r between src and dest nodes.
bool Graph_EdgeExists(const Graph *g, NodeID srcID, NodeID destID, int r) {
	assert(g);
	EdgeID edgeId;
	GrB_Matrix M = Graph_GetRelationMatrix(g, r);
	GrB_Info res = GrB_Matrix_extractElement_UINT64(&edgeId, M, destID, srcID);
	return res == GrB_SUCCESS;
}

static inline Entity *_Graph_GetEntity(const DataBlock *entities, EntityID id) {
	return DataBlock_GetItem(entities, id);
}

/* ============= Matrix synchronization and resizing functions =============== */

/* Resize given matrix, such that its number of row and columns
 * matches the number of nodes in the graph. Also, synchronize
 * matrix to execute any pending operations. */
void _MatrixSynchronize(const Graph *g, RG_Matrix rg_matrix) {
	GrB_Matrix m = RG_Matrix_Get_GrB_Matrix(rg_matrix);
	GrB_Index n_rows;
	GrB_Index n_cols;
	GrB_Matrix_nrows(&n_rows, m);
	GrB_Matrix_ncols(&n_cols, m);
	GrB_Index dims = Graph_RequiredMatrixDim(g);

	// If the graph belongs to one thread, we don't need to lock the mutex.
	if(g->_writelocked) {
		if((n_rows != dims) || (n_cols != dims)) {
			assert(GxB_Matrix_resize(m, dims, dims) == GrB_SUCCESS);
		}

		// Writer under write lock, no need to flush pending changes.
		return;
	}
	// Lock the matrix.
	RG_Matrix_Lock(rg_matrix);

	GrB_wait(); // Apply all pending operations.

	// If the matrix has pending operations or requires
	// a resize, enter critical section.
	if((n_rows != dims) || (n_cols != dims)) {
		// Double-check if resize is necessary.
		GrB_Matrix_nrows(&n_rows, m);
		GrB_Matrix_ncols(&n_cols, m);
		dims = Graph_RequiredMatrixDim(g);
		if((n_rows != dims) || (n_cols != dims)) {
			assert(GxB_Matrix_resize(m, dims, dims) == GrB_SUCCESS);
		}
		// Flush changes to matrix.
		_Graph_ApplyPending(m);
	}
	// Unlock matrix mutex.
	_RG_Matrix_Unlock(rg_matrix);
}

/* Resize matrix to node capacity. */
void _MatrixResizeToCapacity(const Graph *g, RG_Matrix matrix) {
	GrB_Matrix m = RG_Matrix_Get_GrB_Matrix(matrix);
	GrB_Index nrows;
	GrB_Index ncols;
	GrB_Matrix_ncols(&ncols, m);
	GrB_Matrix_nrows(&nrows, m);
	GrB_Index cap = _Graph_NodeCap(g);

	// This policy should only be used in a thread-safe context, so no locking is required.
	if(ncols != cap || nrows != cap) {
		assert(GxB_Matrix_resize(m, cap, cap) == GrB_SUCCESS);
	}
}

/* Do not update matrices. */
void _MatrixNOP(const Graph *g, RG_Matrix matrix) {
	return;
}

/* Define the current behavior for matrix creations and retrievals on this graph. */
void Graph_SetMatrixPolicy(Graph *g, MATRIX_POLICY policy) {
	switch(policy) {
	case SYNC_AND_MINIMIZE_SPACE:
		// Default behavior; forces execution of pending GraphBLAS operations
		// when appropriate and sizes matrices to the current node count.
		g->SynchronizeMatrix = _MatrixSynchronize;
		break;
	case RESIZE_TO_CAPACITY:
		// Bulk insertion and creation behavior; does not force pending operations
		// and resizes matrices to the graph's current node capacity.
		g->SynchronizeMatrix = _MatrixResizeToCapacity;
		break;
	case DISABLED:
		// Used when deleting or freeing a graph; forces no matrix updates or resizes.
		g->SynchronizeMatrix = _MatrixNOP;
		break;
	default:
		assert(false);
	}
}

/* Synchronize and resize all matrices in graph. */
void Graph_ApplyAllPending(Graph *g) {
	RG_Matrix M;

	for(int i = 0; i < array_len(g->labels); i ++) {
		M = g->labels[i];
		g->SynchronizeMatrix(g, M);
	}

	for(int i = 0; i < array_len(g->relations); i ++) {
		M = g->relations[i];
		g->SynchronizeMatrix(g, M);
	}

	if(Config_MaintainTranspose()) {
		for(int i = 0; i < array_len(g->t_relations); i ++) {
			M = g->t_relations[i];
			g->SynchronizeMatrix(g, M);
		}
	}
}

/* ================================ Graph API ================================ */
Graph *Graph_New(size_t node_cap, size_t edge_cap) {
	node_cap = MAX(node_cap, GRAPH_DEFAULT_NODE_CAP);
	edge_cap = MAX(node_cap, GRAPH_DEFAULT_EDGE_CAP);

	Graph *g = rm_malloc(sizeof(Graph));
	g->nodes = DataBlock_New(node_cap, sizeof(Entity), (fpDestructor)FreeEntity);
	g->edges = DataBlock_New(edge_cap, sizeof(Entity), (fpDestructor)FreeEntity);
	g->labels = array_new(RG_Matrix, GRAPH_DEFAULT_LABEL_CAP);
	g->relations = array_new(RG_Matrix, GRAPH_DEFAULT_RELATION_TYPE_CAP);
	g->adjacency_matrix = RG_Matrix_New(GrB_BOOL, node_cap, node_cap);
	g->_t_adjacency_matrix = RG_Matrix_New(GrB_BOOL, node_cap, node_cap);
	g->_zero_matrix = RG_Matrix_New(GrB_BOOL, node_cap, node_cap);
	// If we're maintaining transposed relation matrices, allocate a new array, otherwise NULL-set the pointer.
	g->t_relations = Config_MaintainTranspose() ?
					 array_new(RG_Matrix, GRAPH_DEFAULT_RELATION_TYPE_CAP) : NULL;

	// Initialize a read-write lock scoped to the individual graph
	assert(pthread_rwlock_init(&g->_rwlock, NULL) == 0);
	g->_writelocked = false;

	// Force GraphBLAS updates and resize matrices to node count by default
	Graph_SetMatrixPolicy(g, SYNC_AND_MINIMIZE_SPACE);

	// Synchronization objects initialization.
	assert(pthread_mutex_init(&g->_writers_mutex, NULL) == 0);

	return g;
}

// All graph matrices are required to be squared NXN
// where N is Graph_RequiredMatrixDim.
size_t Graph_RequiredMatrixDim(const Graph *g) {
	// Matrix dimensions should be at least:
	// Number of nodes + number of deleted nodes.
	return g->nodes->itemCount + array_len(g->nodes->deletedIdx);
}

size_t Graph_NodeCount(const Graph *g) {
	assert(g);
	return g->nodes->itemCount;
}

uint Graph_DeletedNodeCount(const Graph *g) {
	assert(g);
	return DataBlock_DeletedItemsCount(g->nodes);
}

size_t Graph_LabeledNodeCount(const Graph *g, int label) {
	GrB_Index nvals = 0;
	GrB_Matrix m = Graph_GetLabelMatrix(g, label);
	if(m) GrB_Matrix_nvals(&nvals, m);
	return nvals;
}

size_t Graph_EdgeCount(const Graph *g) {
	assert(g);
	return g->edges->itemCount;
}

uint Graph_DeletedEdgeCount(const Graph *g) {
	assert(g);
	return DataBlock_DeletedItemsCount(g->edges);
}

int Graph_RelationTypeCount(const Graph *g) {
	return array_len(g->relations);
}

int Graph_LabelTypeCount(const Graph *g) {
	return array_len(g->labels);
}

void Graph_AllocateNodes(Graph *g, size_t n) {
	assert(g);
	DataBlock_Accommodate(g->nodes, n);
}

void Graph_AllocateEdges(Graph *g, size_t n) {
	assert(g);
	DataBlock_Accommodate(g->edges, n);
}

int Graph_GetNode(const Graph *g, NodeID id, Node *n) {
	assert(g);
	n->entity = _Graph_GetEntity(g->nodes, id);
	n->id = id;
	return (n->entity != NULL);
}

int Graph_GetEdge(const Graph *g, EdgeID id, Edge *e) {
	assert(g && id < _Graph_EdgeCap(g));
	e->entity = _Graph_GetEntity(g->edges, id);
	e->id = id;
	return (e->entity != NULL);
}

int Graph_GetNodeLabel(const Graph *g, NodeID nodeID) {
	assert(g);
	int label = GRAPH_NO_LABEL;
	for(int i = 0; i < array_len(g->labels); i++) {
		bool x = false;
		GrB_Matrix M = Graph_GetLabelMatrix(g, i);
		GrB_Info res = GrB_Matrix_extractElement_BOOL(&x, M, nodeID, nodeID);
		if(res == GrB_SUCCESS && x == true) {
			label = i;
			break;
		}
	}

	return label;
}

int Graph_GetEdgeRelation(const Graph *g, Edge *e) {
	assert(g && e);
	NodeID srcNodeID = Edge_GetSrcNodeID(e);
	NodeID destNodeID = Edge_GetDestNodeID(e);
	EdgeID id = ENTITY_GET_ID(e);

	// Search for relation mapping matrix M, where
	// M[dest,src] == edge ID.
	uint relationship_count = array_len(g->relations);
	for(uint i = 0; i < relationship_count; i++) {
		EdgeID edgeId = 0;
		GrB_Matrix M = Graph_GetRelationMatrix(g, i);
		GrB_Info res = GrB_Matrix_extractElement_UINT64(&edgeId, M, srcNodeID, destNodeID);
		if(res != GrB_SUCCESS) continue;

		if(SINGLE_EDGE(edgeId)) {
			EdgeID curEdgeID = SINGLE_EDGE_ID(edgeId);
			if(curEdgeID == id) {
				Edge_SetRelationID(e, i);
				return i;
			}
		} else {
			/* Multiple edges exists between src and dest
			 * see if given edge is one of them. */
			EdgeID *edges = (EdgeID *)edgeId;
			int edge_count = array_len(edges);
			for(int j = 0; j < edge_count; j++) {
				if(edges[j] == id) {
					Edge_SetRelationID(e, i);
					return i;
				}
			}
		}
	}

	// We must be able to find edge relation.
	assert(false);
	return GRAPH_NO_RELATION;
}

void Graph_GetEdgesConnectingNodes(const Graph *g, NodeID srcID, NodeID destID, int r,
								   Edge **edges) {
	assert(g && r < Graph_RelationTypeCount(g) && edges);

	// Invalid relation type specified; this can occur on multi-type traversals like:
	// MATCH ()-[:real_type|fake_type]->()
	if(r == GRAPH_UNKNOWN_RELATION) return;

	Node srcNode = GE_NEW_NODE();
	Node destNode = GE_NEW_NODE();
	assert(Graph_GetNode(g, srcID, &srcNode));
	assert(Graph_GetNode(g, destID, &destNode));

	if(r != GRAPH_NO_RELATION) {
		_Graph_GetEdgesConnectingNodes(g, srcID, destID, r, edges);
	} else {
		// Relation type missing, scan through each edge type.
		int relationCount = Graph_RelationTypeCount(g);
		for(int i = 0; i < relationCount; i++) {
			_Graph_GetEdgesConnectingNodes(g, srcID, destID, i, edges);
		}
	}
}

void Graph_CreateNode(Graph *g, int label, Node *n) {
	assert(g);

	NodeID id;
	Entity *en = DataBlock_AllocateItem(g->nodes, &id);
	n->id = id;
	n->entity = en;
	en->prop_count = 0;
	en->properties = NULL;

	if(label != GRAPH_NO_LABEL) {
		// Try to set matrix at position [id, id]
		// incase of a failure, scale matrix.
		RG_Matrix matrix = g->labels[label];
		GrB_Matrix m = RG_Matrix_Get_GrB_Matrix(matrix);
		GrB_Info res = GrB_Matrix_setElement_BOOL(m, true, id, id);
		if(res != GrB_SUCCESS) {
			_MatrixResizeToCapacity(g, matrix);
			assert(GrB_Matrix_setElement_BOOL(m, true, id, id) == GrB_SUCCESS);
		}
	}
}

static void add_relation_id (EdgeID* z_out, EdgeID x_in, EdgeID y_newid)
{
	EdgeID *ids;
	/* Single edge ID,
	 * switching from single edge ID to multiple IDs. */
    printf("    add_relation_id: add edgeID 0x%x\n", y_newid ); fflush(stdout);
	if(SINGLE_EDGE(x_in)) {
        printf("        add_relation_id: before: single edge\n", z_out, array_len(z_out)); fflush(stdout);
		ids = array_new(EdgeID, 2);
		ids = array_append(ids, SINGLE_EDGE_ID(x_in));
		ids = array_append(ids, SINGLE_EDGE_ID(y_newid));
		// TODO: Make sure MSB of ids isn't on.
		*z_out = (EdgeID)ids;
	} else {
        printf("        add_relation_id: before: edgeID %p (%ld), count %ld\n", x_in, (long) x_in, array_len(x_in)); fflush(stdout);
		// Multiple edges, adding another edge.
		ids = (EdgeID *)(x_in);
		ids = array_append(ids, SINGLE_EDGE_ID(y_newid));
		*z_out = (EdgeID)ids;
	}
    printf("        add_relation_id: after: add edgeID %p (%ld), count %ld\n", *z_out, (long) *z_out, array_len(*z_out)); fflush(stdout);
}

void Graph_FormConnection(Graph *g, NodeID src, NodeID dest, EdgeID edge_id, int r) {
	GrB_Matrix adj = Graph_GetAdjacencyMatrix(g);
	GrB_Matrix tadj = Graph_GetTransposedAdjacencyMatrix(g);
	GrB_Matrix relationMat = Graph_GetRelationMatrix(g, r);
        GrB_Info info;

        printf("FORM CONNECTION: src %ld, dest %ld, edge id %ld, r %u\n",
                src, dest, edge_id, r); fflush(stdout);

	// Rows represent source nodes, columns represent destination nodes.
	GrB_Matrix_setElement_BOOL(adj, true, src, dest);
	GrB_Matrix_setElement_BOOL(tadj, true, dest, src);
	edge_id = SET_MSB(edge_id);

        EdgeID relationMat_ij;
        info = GrB_Matrix_extractElement_UINT64 ((uint64_t*)&relationMat_ij, relationMat, src, dest);
        if (info == GrB_NO_VALUE) {
            printf("  Add new edge_id %p (%lu)", edge_id, (long) edge_id); fflush(stdout);
             relationMat_ij = edge_id;
        } else {
             INFO_SUCCESS(info);
             add_relation_id (&relationMat_ij, relationMat_ij, edge_id);
        }
        info = GrB_Matrix_setElement_UINT64 (relationMat, (uint64_t)relationMat_ij, src, dest);
    GxB_Matrix_fprint(relationMat, "relationMat", (GxB_Print_Level) 3, stdout); fflush(stdout);
	INFO_SUCCESS(info);

	// Update the transposed matrix if one is present.
	if(Config_MaintainTranspose()) {

		// Perform the same update to the J,I coordinates of the transposed matrix.
		GrB_Matrix t_relationMat = Graph_GetTransposedRelationMatrix(g, r);
                EdgeID relationMat_ji;
                info = GrB_Matrix_extractElement_UINT64 ((uint64_t*)&relationMat_ji, t_relationMat, dest, src);
                if (info == GrB_NO_VALUE) {
                    printf("  TRANSPOSE: add new edge_id %p (%lu)",
                            edge_id, (long) edge_id); fflush(stdout);
                     relationMat_ji = edge_id;
                } else {
                     INFO_SUCCESS(info);
                     add_relation_id (&relationMat_ji, relationMat_ji, edge_id);
                }
                info = GrB_Matrix_setElement_UINT64 (t_relationMat, (uint64_t)relationMat_ji, dest, src);
                INFO_SUCCESS(info);
        GxB_Matrix_fprint(t_relationMat, "t_relationMat", (GxB_Print_Level) 3, stdout); fflush(stdout);
	}
}

int Graph_ConnectNodes(Graph *g, NodeID src, NodeID dest, int r, Edge *e) {
	Node srcNode = GE_NEW_NODE();
	Node destNode = GE_NEW_NODE();

	assert(Graph_GetNode(g, src, &srcNode));
	assert(Graph_GetNode(g, dest, &destNode));
	assert(g && r < Graph_RelationTypeCount(g));

	EdgeID id;
	Entity *en = DataBlock_AllocateItem(g->edges, &id);
	en->prop_count = 0;
	en->properties = NULL;
	e->id = id;
	e->entity = en;
	e->relationID = r;
	e->srcNodeID = src;
	e->destNodeID = dest;
	Graph_FormConnection(g, src, dest, id, r);
	return 1;
}

/* Retrieves all either incoming or outgoing edges
 * to/from given node N, depending on given direction. */
void Graph_GetNodeEdges(const Graph *g, const Node *n, GRAPH_EDGE_DIR dir, int edgeType,
						Edge **edges) {
	assert(g && n && edges);
        GrB_Descriptor desc;
	GrB_Matrix M;
        GrB_Vector v;
	NodeID srcNodeID;
	NodeID destNodeID;
        GrB_Index nrows, ncols;
        GrB_Info info;

	if(edgeType == GRAPH_UNKNOWN_RELATION) return;

        /* If a relationship type is specified, retrieve the appropriate relation matrix;
         * otherwise use the overall adjacency matrix. */
        if(edgeType == GRAPH_NO_RELATION) M = Graph_GetAdjacencyMatrix(g);
        else M = Graph_GetRelationMatrix(g, edgeType);

        GrB_Descriptor_new (&desc) ;
        GrB_Descriptor_set (desc, GrB_OUTP, GrB_REPLACE) ;

        // Outgoing.
	if(dir == GRAPH_EDGE_DIR_OUTGOING || dir == GRAPH_EDGE_DIR_BOTH) {
		srcNodeID = ENTITY_GET_ID(n);
                GrB_Matrix_nrows (&nrows, M);
                GrB_Vector_new (&v, GrB_BOOL, nrows);
                GrB_Vector_setElement_BOOL (v, true, srcNodeID);
                GrB_vxm (v, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL, v, M, desc) ;

                GrB_Type vtype;
                GrB_Index n, nvals, *I;
                bool *X;
                GrB_wait();
                GxB_Vector_export (&v, &vtype, &n, &nvals, &I, (void**)&X, GrB_NULL);
                assert (vtype == GrB_BOOL);

                for (GrB_Index k = 0; k < nvals; ++k)
                        Graph_GetEdgesConnectingNodes(g, srcNodeID, I[k], edgeType, edges);

                free (X);
                free (I);
                // v is freed in GxB_Vector_export
	}

	// Incoming.
	if(dir == GRAPH_EDGE_DIR_INCOMING || dir == GRAPH_EDGE_DIR_BOTH) {
                GrB_Descriptor_set (desc, GrB_INP1, GrB_TRAN);
		destNodeID = ENTITY_GET_ID(n);
                GrB_Matrix_ncols (&ncols, M);
                GrB_Vector_new (&v, GrB_BOOL, ncols);
                GrB_Vector_setElement_BOOL (v, true, destNodeID);
                GrB_vxm (v, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL, v, M, desc) ;

                GrB_Type vtype;
                GrB_Index n, nvals, *I;
                bool *X;
                GrB_wait();
                GxB_Vector_export (&v, &vtype, &n, &nvals, &I, (void**)&X, GrB_NULL);
                assert (vtype == GrB_BOOL);

                for (GrB_Index k = 0; k < nvals; ++k)
                        Graph_GetEdgesConnectingNodes(g, I[k], destNodeID, edgeType, edges);

                free (X);
                free (I);
                // v is freed in GxB_Vector_export
	}

        GrB_free (&desc);
}

/* Removes an edge from Graph and updates graph relevent matrices. */
int Graph_DeleteEdge(Graph *g, Edge *e) {
	uint64_t x;
	GrB_Matrix R;
	GrB_Matrix M;
	GrB_Info info;
	EdgeID edge_id;
	GrB_Matrix TR = GrB_NULL;
	int r = Edge_GetRelationID(e);
	NodeID src_id = Edge_GetSrcNodeID(e);
	NodeID dest_id = Edge_GetDestNodeID(e);

	R = Graph_GetRelationMatrix(g, r);
	if(Config_MaintainTranspose()) TR = Graph_GetTransposedRelationMatrix(g, r);

	// Test to see if edge exists.
	info = GrB_Matrix_extractElement(&edge_id, R, src_id, dest_id);
	if(info != GrB_SUCCESS) return 0;

	if(SINGLE_EDGE(edge_id)) {
		// Single edge of type R connecting src to dest, delete entry.
		assert(GxB_Matrix_Delete(R, src_id, dest_id) == GrB_SUCCESS);
		if(TR) assert(GxB_Matrix_Delete(TR, dest_id, src_id) == GrB_SUCCESS);

		// See if source is connected to destination with additional edges.
		bool connected = false;
		int relationCount = Graph_RelationTypeCount(g);
		for(int i = 0; i < relationCount; i++) {
			if(i == r) continue;
			M = Graph_GetRelationMatrix(g, i);
			info = GrB_Matrix_extractElement(&x, M, src_id, dest_id);
			if(info == GrB_SUCCESS) {
				connected = true;
				break;
			}
		}

		/* There are no additional edges connecting source to destination
		 * Remove edge from THE adjacency matrix. */
		if(!connected) {
			M = Graph_GetAdjacencyMatrix(g);
			assert(GxB_Matrix_Delete(M, src_id, dest_id) == GrB_SUCCESS);

			M = Graph_GetTransposedAdjacencyMatrix(g);
			assert(GxB_Matrix_Delete(M, dest_id, src_id) == GrB_SUCCESS);
		}
	} else {
		/* Multiple edges connecting src to dest
		 * locate specific edge and remove it
		 * revert back from array representation to edge ID
		 * incase we're left with a single edge connecting src to dest. */

		int i = 0;
		EdgeID id = ENTITY_GET_ID(e);
		EdgeID *edges = (EdgeID *)edge_id;
		int edge_count = array_len(edges);

		// Locate edge within edge array.
		for(; i < edge_count; i++) if(edges[i] == id) break;
		assert(i < edge_count);

		/* Remove edge from edge array
		 * migrate last edge ID and reduce array size.
		 * TODO: reallocate array of size / capacity ratio is high. */
		edges[i] = edges[edge_count - 1];
		array_pop(edges);

		/* Incase we're left with a single edge connecting src to dest
		 * revert back from array to scalar. */
		if(array_len(edges) == 1) {
			edge_id = edges[0];
			array_free(edges);
			GrB_Matrix_setElement(R, SET_MSB(edge_id), src_id, dest_id);
		}

		if(TR) {
			/* We must make the matching updates to the transposed matrix.
			 * First, extract the element that is known to be an edge array. */
			info = GrB_Matrix_extractElement(edges, TR, dest_id, src_id);
			INFO_SUCCESS(info);
			// Replace the deleted edge with the last edge in the matrix.
			edges[i] = edges[edge_count - 1];
			array_pop(edges);
			// Free and replace the array if it now has 1 element.
			if(array_len(edges) == 1) {
				edge_id = edges[0];
				array_free(edges);
				GrB_Matrix_setElement(TR, SET_MSB(edge_id), dest_id, src_id);
			}
		}
	}

	// Free and remove edges from datablock.
	DataBlock_DeleteItem(g->edges, ENTITY_GET_ID(e));
	return 1;
}

void Graph_DeleteNode(Graph *g, Node *n) {
	/* Assumption, node is completely detected,
	 * there are no incoming nor outgoing edges
	 * leading to / from node. */
	assert(g && n);

	// Clear label matrix at position node ID.
	uint32_t label_count = array_len(g->labels);
	for(int i = 0; i < label_count; i++) {
		GrB_Matrix M = Graph_GetLabelMatrix(g, i);
		GxB_Matrix_Delete(M, ENTITY_GET_ID(n), ENTITY_GET_ID(n));
	}

	DataBlock_DeleteItem(g->nodes, ENTITY_GET_ID(n));
}

static void free_edge (Graph *g, EdgeID id)
{
     if (SINGLE_EDGE(id)) {
         printf("    free_edge: free single edge %ld\n", SINGLE_EDGE_ID(id)); fflush(stdout);
          DataBlock_DeleteItem (g->edges, SINGLE_EDGE_ID(id));
     } else {

         EdgeID *ids = (EdgeID *) id;
         if (*ids == NULL)
             return;
         uint id_count = array_len(ids);
         printf("    free_edge: free multi edge %p, nedges 0x%x\n", ids, id_count);
         fflush(stdout);
         for (uint i = 0; i < id_count; i++) {
             DataBlock_DeleteItem(g->edges, ids[i]);
         }
         if (id_count > 0) {
             array_free(ids);
             *ids = NULL;
         } else {
             *ids = NULL;
         }
         printf("    free_edge: After free: ids %p, *ids %p\n", ids, *ids); fflush(stdout);
     }
     //
}

static void _Graph_FreeRelationMatrixEdges(Graph *g, GrB_Matrix grb_matrix, GrB_Index *nvals_alloc_ptr, EdgeID **valptr)
{
        GrB_Index nvals_alloc = *nvals_alloc_ptr;
        EdgeID *val = *valptr;

        GrB_Index nvals, nvals_saved;
        GrB_Info info;

        info = GrB_Matrix_nvals (&nvals, grb_matrix);
        INFO_SUCCESS(info);

        if (nvals == 0) {
            printf("   nvals == 0, empty relationship matrix\n"); fflush(stdout);
            return;
        }

    GxB_Matrix_fprint(grb_matrix, "grb_matrix", (GxB_Print_Level) 3, stdout);
    // why are we allocating new edges?
        printf(" Nvals alloc %ld, nvals %ld\n", nvals_alloc, nvals); fflush(stdout);
        if (nvals > nvals_alloc) {
             EdgeID *tmp_e;
             tmp_e = realloc (val, nvals * sizeof (*val));
             assert(tmp_e);
             val = tmp_e;
        }
        printf("  After if nvals > nvals_alloc\n"); fflush(stdout);

#if !defined(NDEBUG)
        nvals_saved = nvals;
#endif

        // Extract tuples into val from grb_matrix
        info = GrB_Matrix_extractTuples (GrB_NULL, GrB_NULL, (uintptr_t*)val, &nvals, grb_matrix);
        INFO_SUCCESS(info);
#if !defined(NDEBUG)
        assert (nvals == nvals_saved);
#endif

    for (GrB_Index k = 0; k < nvals; ++k) {
        printf("val[%ld] %p\n", k, val[k]); fflush(stdout);
    }

        //OMP(omp parallel for)
        // Issue appears to be here. Why?
    for (GrB_Index k = 0; k < nvals; ++k)
                  free_edge (g, val[k]);
}

static void _Graph_FreeRelationMatrices(Graph *g) {
	uint relationCount = Graph_RelationTypeCount(g);
        GrB_Index nvals_alloc = 0;
        EdgeID *val = NULL;
	for(uint i = 0; i < relationCount; i++) {
		RG_Matrix M = g->relations[i];

                GrB_Matrix grb_matrix = M->grb_matrix;
                _Graph_FreeRelationMatrixEdges (g, grb_matrix, &nvals_alloc, &val);

		// Free the matrix itself.
		RG_Matrix_Free(M);

		// Perform the same update to transposed matrices.
		if(Config_MaintainTranspose()) {
			RG_Matrix TM = g->t_relations[i];

                        grb_matrix = TM->grb_matrix;
                        _Graph_FreeRelationMatrixEdges (g, grb_matrix, &nvals_alloc, &val);

			// Free the matrix itself.
			RG_Matrix_Free(TM);
		}
	}
        free (val);
}

static void _BulkDeleteNodes(Graph *g, Node *nodes, uint node_count,
							 uint *node_deleted, uint *edge_deleted) {
	assert(g && g->_writelocked && nodes && node_count > 0);

	/* Create a matrix M where M[j,i] = 1 where:
	 * Node i is connected to node j. */

    GrB_Info info;
    /* XXX: Need to count number of edges deleted for some reason */
    uint64_t nedge_deleted, nedge_deleted_row;
    GrB_Matrix tmp;
    GrB_Matrix tmp_transpose;
    GrB_Index dim;

	GrB_Matrix adj;                     // Adjacency matrix.
	GrB_Matrix tadj;                    // Transposed adjacency matrix.

	//GrB_Index *node_idx = NULL;
    GrB_Index *node_idx = (GrB_Index *) malloc(node_count * sizeof(GrB_Index));

	adj = Graph_GetAdjacencyMatrix(g);
	tadj = Graph_GetTransposedAdjacencyMatrix(g);
#if 1
    GrB_Index nrows;
    GrB_Matrix_nrows(&nrows, adj);
    GrB_Index ncols;
    GrB_Matrix_ncols(&ncols, adj);
    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, adj);

    printf(" On entry: edge_deleted %d, node_deleted %d\n", *edge_deleted, *node_deleted); fflush(stdout);
    for (long i=0; i<node_count; i++) {
        printf(" nodes[%ld] = %p\n", i, nodes+i); fflush(stdout);
    }

    printf("  Initial values for adj and tadj\n");
    GxB_Matrix_fprint(adj, "adj", (GxB_Print_Level) 3, stdout);
    GxB_Matrix_fprint(tadj, "tadj", (GxB_Print_Level) 3, stdout);
#endif

    // Get the node ids
    dim = Graph_RequiredMatrixDim (g);
    printf("  adj mat: dim %ld, nrows %ld ncols %ld nvals %ld\n", dim, nrows, ncols, nvals); fflush(stdout);

    for(uint i = 0; i < node_count; i++) {
        //node_idx[i] = ENTITY_GET_ID(nodes + i);
        node_idx[i] = ENTITY_GET_ID(&nodes[i]);
        printf("  node_idx[%d] = %ld\n", i, node_idx[i]);
    }

    // Create tmp of size dim x node_count
    info = GrB_Matrix_new (&tmp, GrB_UINT64, dim, node_count);
    INFO_SUCCESS(info);

    // Extract node_idx values from adj matrix into tmp
    // tmp will hold the values that need to be deleted from adj matrix
    // Note that it is based on I, J so tmp(0,3) means A(I(0), J(3)) = A(0, 1)
    info = GrB_extract (tmp, GrB_NULL, GrB_NULL, adj, GrB_ALL, dim, node_idx, node_count,
                            GrB_NULL);
    INFO_SUCCESS(info);
    GrB_Matrix_nvals (&nedge_deleted, tmp);
    printf(" tmp nvals %ld\n", nedge_deleted);
    GxB_Matrix_fprint(tmp, "tmp", (GxB_Print_Level) 3, stdout);

    // Create tmp transpose
    info = GrB_Matrix_new (&tmp_transpose, GrB_UINT64, node_count, dim);
    INFO_SUCCESS(info);
    info = GrB_transpose(tmp_transpose, NULL, NULL, tmp, NULL);
    INFO_SUCCESS(info);
    GxB_Matrix_fprint(tmp_transpose, "tmp_transpose", (GxB_Print_Level) 3, stdout);

    // Assign tmp back into adj so that these (now zero) values will be deleted
    GrB_Matrix_clear (tmp); // Zero out the values
    info = GrB_assign (adj, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, node_idx, node_count,
                           GrB_NULL);
    INFO_SUCCESS(info);
    printf("  After assign tmp back into adj\n");
    GxB_Matrix_fprint(adj, "adj", (GxB_Print_Level) 3, stdout);

    // Assign tmp back into tadj so that these (now zero) values will be deleted
    GrB_Matrix_clear (tmp_transpose); // Zero out the values
    //info = GrB_assign (tadj, GrB_NULL, GrB_NULL, tmp_transpose, GrB_ALL, dim, node_idx, node_count,
    info = GrB_assign (tadj, GrB_NULL, GrB_NULL, tmp_transpose, node_idx, node_count, GrB_ALL, dim,
                       GrB_NULL);
    INFO_SUCCESS(info);
    printf("  Did tadj assign actually work?\n"); fflush(stdout);
    GxB_Matrix_fprint(tadj, "tadj", (GxB_Print_Level) 3, stdout);


#if 0  // Rewrite the transpose update above
    //info = GxB_Matrix_resize (tmp, dim, node_count);
    //INFO_SUCCESS(info);
    info = GrB_extract (tmp, GrB_NULL, GrB_NULL, adj, GrB_ALL, dim, node_idx, node_count,
                        GrB_NULL);
    INFO_SUCCESS(info);
    GrB_Matrix_nvals (&nedge_deleted_row, tmp);
    nedge_deleted += nedge_deleted_row;

    GrB_Matrix_clear (tmp);
    info = GrB_Matrix_assign (adj, GrB_NULL, GrB_NULL, tmp, node_idx, node_count, GrB_ALL, dim,
                              GrB_NULL);
    INFO_SUCCESS(info);


    info = GrB_Matrix_assign (tadj, GrB_NULL, GrB_NULL, tmp, node_idx, node_count, GrB_ALL, dim,
                              GrB_NULL);
    INFO_SUCCESS(info);
#endif

    // FIXME: This computation isn't right. Should be 2 nodes and 4 edges deleted.
	*node_deleted += node_count; // XXX: Um... What if it already was cleared?
	*edge_deleted += nedge_deleted;

	printf("  Before relationship matrices, edge_deleted %d, node deleted %d\n",
	        *edge_deleted, *node_deleted); fflush(stdout);

	//////////////////////////////////////////////////////////////////////
    // Free and remove implicit edges from relation matrices.
	GrB_Index nvals_alloc = 0;
	EdgeID *val = NULL;
	const int relation_count = Graph_RelationTypeCount(g);
	const bool update_transpose = Config_MaintainTranspose();
	for(int i = 0; i < relation_count; i++) {
	    printf(" In loop for relation count: i %d\n", i); fflush(stdout);
		GrB_Matrix R = Graph_GetRelationMatrix(g, i);


        info = GxB_Matrix_resize (tmp, dim, node_count);
        INFO_SUCCESS(info);

        info = GrB_extract (tmp, GrB_NULL, GrB_NULL, R, GrB_ALL, dim, node_idx, node_count,
                            GrB_NULL);
        INFO_SUCCESS(info);

        printf("  Before Free Relation matrix edges\n"); fflush(stdout);
        //_Graph_FreeRelationMatrixEdges (g, tmp, &nvals_alloc, &val);
        _Graph_FreeRelationMatrixEdges (g, tmp, &nvals_alloc, &val);
        printf("  After Free Relation matrix edges\n"); fflush(stdout);

        GrB_Matrix Rt = GrB_NULL;
        info = GxB_Matrix_resize (tmp_transpose, node_count, dim);
        if (update_transpose) {
            Rt = Graph_GetTransposedRelationMatrix(g, i);
            // Issue here - I think tmp is NULL
            info = GrB_extract (tmp_transpose, GrB_NULL, GrB_NULL, Rt, node_idx, node_count,
                                 GrB_ALL, dim, GrB_NULL);
            INFO_SUCCESS(info);
            // XXX: Apparently the transpose may hold different edge IDs?
            _Graph_FreeRelationMatrixEdges (g, tmp_transpose, &nvals_alloc, &val);
        }
        printf("  After update transpose\n"); fflush(stdout);

        //// WTH is this stuff? We already free'd tmp
        ////
        GrB_Matrix_clear (tmp);
        info = GrB_assign (R, GrB_NULL, GrB_NULL, tmp, GrB_ALL, dim, node_idx, node_count,
                           GrB_NULL);
        INFO_SUCCESS(info);
        if (update_transpose) {
             info = GrB_assign (Rt, GrB_NULL, GrB_NULL, tmp_transpose, node_idx, node_count,
                     GrB_ALL, dim, GrB_NULL);
             INFO_SUCCESS(info);
        }

        info = GxB_Matrix_resize (tmp, node_count, dim);
        INFO_SUCCESS(info);
        info = GrB_extract (tmp, GrB_NULL, GrB_NULL, R, node_idx, node_count,
                GrB_ALL, dim, GrB_NULL);
        INFO_SUCCESS(info);

        _Graph_FreeRelationMatrixEdges (g, tmp, &nvals_alloc, &val);
        if (update_transpose) {
            info = GxB_Matrix_resize (tmp_transpose, dim, node_count);
            INFO_SUCCESS(info);
             info = GrB_extract (tmp_transpose, GrB_NULL, GrB_NULL, Rt,
                                 GrB_ALL, dim, node_idx, node_count, GrB_NULL);
             INFO_SUCCESS(info);
             // XXX: Apparently the transpose may hold different edge IDs?
             _Graph_FreeRelationMatrixEdges (g, tmp, &nvals_alloc, &val);
        }
        ////
        printf(" After free edges 2\n"); fflush(stdout);

        GrB_Matrix_clear (tmp);
        info = GrB_assign (R, GrB_NULL, GrB_NULL, tmp, node_idx, node_count,  GrB_ALL, dim,
                           GrB_NULL);
        INFO_SUCCESS(info);
        if (update_transpose) {
             info = GrB_assign (Rt, GrB_NULL, GrB_NULL, tmp_transpose, GrB_ALL, dim, node_idx, node_count,
                                GrB_NULL);
             INFO_SUCCESS(info);
        }
        printf(" After assign 2\n"); fflush(stdout);
	}

        info = GxB_Matrix_resize (tmp, node_count, node_count);
        INFO_SUCCESS(info);

        printf(" Before delete nodes\n"); fflush(stdout);
	/* Delete nodes
	 * All nodes marked for deleteion are detected, no incoming / outgoing edges. */
	int node_type_count = Graph_LabelTypeCount(g);
	for(int i = 0; i < node_type_count; i++) {
                GrB_Matrix L = Graph_GetLabelMatrix(g, i); // XXX: Apparently all diagonal.
		GrB_Matrix_assign (L, GrB_NULL, GrB_NULL, tmp, node_idx, node_count, node_idx, node_count,
                                   GrB_NULL);
	}

	for(uint i = 0; i < node_count; i++) {
		Node *n = nodes + i;
		DataBlock_DeleteItem(g->nodes, ENTITY_GET_ID(n));
	}

	// Clean up.
        free (val);
	GrB_free (&tmp);
}

static void _BulkDeleteEdges(Graph *g, Edge *edges, size_t edge_count) {
	assert(g && g->_writelocked && edges && edge_count > 0);

	int relationCount = Graph_RelationTypeCount(g);
	GrB_Matrix masks[relationCount];
	for(int i = 0; i < relationCount; i++) masks[i] = NULL;
	bool update_adj_matrices = false;

	for(int i = 0; i < edge_count; i++) {
		Edge *e = edges + i;
		int r = Edge_GetRelationID(e);
		NodeID src_id = Edge_GetSrcNodeID(e);
		NodeID dest_id = Edge_GetDestNodeID(e);
		EdgeID edge_id;
		GrB_Matrix R = Graph_GetRelationMatrix(g, r);  // Relation matrix.
		GrB_Matrix TR = Config_MaintainTranspose() ? Graph_GetTransposedRelationMatrix(g, r) : NULL;
		GrB_Matrix_extractElement(&edge_id, R, src_id, dest_id);

		if(SINGLE_EDGE(edge_id)) {
			update_adj_matrices = true;
			GrB_Matrix mask = masks[r];    // mask noteing all deleted edges.
			// Get mask of this relation type.
			if(mask == NULL) {
				GrB_Matrix_new(&mask, GrB_BOOL, Graph_RequiredMatrixDim(g), Graph_RequiredMatrixDim(g));
				masks[r] = mask;
			}
			// Update mask.
			GrB_Matrix_setElement_BOOL(mask, true, src_id, dest_id);
		} else {
			/* Multiple edges connecting src to dest
			 * locate specific edge and remove it
			 * revert back from array representation to edge ID
			 * incase we're left with a single edge connecting src to dest. */

			int i = 0;
			EdgeID id = ENTITY_GET_ID(e);
			EdgeID *multi_edges = (EdgeID *)edge_id;
			int multi_edge_count = array_len(multi_edges);

			// Locate edge within edge array.
			for(; i < multi_edge_count; i++) if(multi_edges[i] == id) break;
			assert(i < multi_edge_count);

			/* Remove edge from edge array
			 * migrate last edge ID and reduce array size.
			 * TODO: reallocate array of size / capacity ratio is high. */
			multi_edges[i] = multi_edges[multi_edge_count - 1];
			array_pop(multi_edges);

			/* Incase we're left with a single edge connecting src to dest
			 * revert back from array to scalar. */
			if(array_len(multi_edges) == 1) {
				edge_id = multi_edges[0];
				array_free(multi_edges);
				GrB_Matrix_setElement(R, SET_MSB(edge_id), src_id, dest_id);
			}

			if(TR) {
				/* We must make the matching updates to the transposed matrix.
				 * First, extract the element that is known to be an edge array. */
				GrB_Matrix_extractElement(&edge_id, TR, dest_id, src_id);
				multi_edges = (EdgeID *)edge_id;
				int multi_edge_count = array_len(multi_edges);
				multi_edges[i] = multi_edges[multi_edge_count - 1];
				array_pop(multi_edges);
				// Free and replace the array if it now has 1 element.
				if(array_len(multi_edges) == 1) {
					edge_id = multi_edges[0];
					array_free(multi_edges);
					GrB_Matrix_setElement(TR, SET_MSB(edge_id), dest_id, src_id);
				}
			}
		}

		// Free and remove edges from datablock.
		DataBlock_DeleteItem(g->edges, ENTITY_GET_ID(e));
	}

	if(update_adj_matrices) {
		GrB_Matrix remaining_mask;
		GrB_Matrix_new(&remaining_mask, GrB_BOOL, Graph_RequiredMatrixDim(g), Graph_RequiredMatrixDim(g));
		GrB_Descriptor desc;    // GraphBLAS descriptor.
		GrB_Descriptor_new(&desc);
		// Descriptor sets to clear entry according to mask.
		GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP);

		// Clear updated output matrix before assignment.
		GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

		for(int r = 0; r < relationCount; r++) {
			GrB_Matrix mask = masks[r];
			GrB_Matrix R = Graph_GetRelationMatrix(g, r);  // Relation matrix.
			if(mask) {
				// Remove every entry of R marked by Mask.
				// Desc: GrB_MASK = GrB_COMP,  GrB_OUTP = GrB_REPLACE.
				// R = R & !mask.
				GrB_Matrix_apply(R, mask, GrB_NULL, GrB_IDENTITY_UINT64, R, desc);
				if(Config_MaintainTranspose()) {
					GrB_Matrix tM = Graph_GetTransposedRelationMatrix(g, r);  // Transposed relation mapping matrix.
					// Transpose mask (this cannot be done by descriptor).
					GrB_transpose(mask, GrB_NULL, GrB_NULL, mask, GrB_NULL);
					// tM = tM & !mask.
					GrB_Matrix_apply(tM, mask, GrB_NULL, GrB_IDENTITY_UINT64, tM, desc);
				}
				GrB_free(&mask);
			}

			// Collect remaining edges. remaining_mask = remaining_mask + R.
			GrB_eWiseAdd_Matrix_Semiring(remaining_mask, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, remaining_mask,
										 R, GrB_NULL);
		}

		GrB_Matrix adj_matrix = Graph_GetAdjacencyMatrix(g);
		GrB_Matrix t_adj_matrix = Graph_GetTransposedAdjacencyMatrix(g);
		// To calculate edges to delete, remove all the remaining edges from "The" adjency matrix.
		// Set descriptor mask to default.
		GrB_Descriptor_set(desc, GrB_MASK, GxB_DEFAULT);
		// adj_matrix = adj_matrix & remaining_mask.
		GrB_Matrix_apply(adj_matrix, remaining_mask, GrB_NULL, GrB_IDENTITY_BOOL, adj_matrix, desc);
		// Transpose remaining_mask.
		GrB_transpose(remaining_mask, GrB_NULL,  GrB_NULL, remaining_mask, GrB_NULL);
		// t_adj_matrix = t_adj_matrix & remaining_mask.
		GrB_Matrix_apply(t_adj_matrix, remaining_mask, GrB_NULL, GrB_IDENTITY_BOOL, t_adj_matrix, desc);

		GrB_free(&remaining_mask);
		GrB_free(&desc);
	}
}

/* Removes both nodes and edges from graph. */
void Graph_BulkDelete(Graph *g, Node *nodes, uint node_count, Edge *edges, uint edge_count,
					  uint *node_deleted, uint *edge_deleted) {
	assert(g);

	*edge_deleted = 0;
	*node_deleted = 0;

	for (long i=0; i<node_count; i++) {
	    printf("nodes[%ld]: %p, id %ld\n", i, nodes+i, ENTITY_GET_ID(nodes+i)); fflush(stdout);
	}

	if(node_count) _BulkDeleteNodes(g, nodes, node_count, node_deleted, edge_deleted);

	if(edge_count) {
		// Filter out explicit edges which were removed by _BulkDeleteNodes.
		if(node_count) {
			for(int i = 0; i < edge_count; i++) {
				Edge *e = edges + i;
				NodeID src = Edge_GetSrcNodeID(e);
				NodeID dest = Edge_GetDestNodeID(e);

				if(!DataBlock_GetItem(g->nodes, src) || !DataBlock_GetItem(g->nodes, dest)) {
					/* Edge already removed due to node removal.
					* Replace current edge with last edge. */
					edges[i] = edges[edge_count - 1];

					// Update indices.
					i--;
					edge_count--;
				}
			}
		}

		/* it might be that edge_count dropped to 0
		 * due to implicit edge deletion. */
		if(edge_count == 0) return;

		// Removing duplicates.
#define is_edge_lt(a, b) (ENTITY_GET_ID((a)) < ENTITY_GET_ID((b)))
		QSORT(Edge, edges, edge_count, is_edge_lt);

		size_t uniqueIdx = 0;
		for(int i = 0; i < edge_count; i++) {
			// As long as current is the same as follows.
			while(i < edge_count - 1 && ENTITY_GET_ID(edges + i) == ENTITY_GET_ID(edges + i + 1)) i++;

			if(uniqueIdx < i) edges[uniqueIdx] = edges[i];
			uniqueIdx++;
		}

		edge_count = uniqueIdx;
		_BulkDeleteEdges(g, edges, edge_count);
	}

	*edge_deleted += edge_count;
}

DataBlockIterator *Graph_ScanNodes(const Graph *g) {
	assert(g);
	return DataBlock_Scan(g->nodes);
}

DataBlockIterator *Graph_ScanEdges(const Graph *g) {
	assert(g);
	return DataBlock_Scan(g->edges);
}

int Graph_AddLabel(Graph *g) {
	assert(g);
	RG_Matrix m = RG_Matrix_New(GrB_BOOL, Graph_RequiredMatrixDim(g), Graph_RequiredMatrixDim(g));
	array_append(g->labels, m);
	return array_len(g->labels) - 1;
}

int Graph_AddRelationType(Graph *g) {
	assert(g);

	size_t dims = Graph_RequiredMatrixDim(g);
	RG_Matrix m = RG_Matrix_New(GrB_UINT64, dims, dims);
	g->relations = array_append(g->relations, m);
	if(Config_MaintainTranspose()) {
		RG_Matrix tm = RG_Matrix_New(GrB_UINT64, dims, dims);
		g->t_relations = array_append(g->t_relations, tm);
	}

	int relationID = Graph_RelationTypeCount(g) - 1;
	return relationID;
}

GrB_Matrix Graph_GetAdjacencyMatrix(const Graph *g) {
	assert(g);
	RG_Matrix m = g->adjacency_matrix;
	g->SynchronizeMatrix(g, m);
	return RG_Matrix_Get_GrB_Matrix(m);
}

// Get the transposed adjacency matrix.
GrB_Matrix Graph_GetTransposedAdjacencyMatrix(const Graph *g) {
	assert(g);
	RG_Matrix m = g->_t_adjacency_matrix;
	g->SynchronizeMatrix(g, m);
	return RG_Matrix_Get_GrB_Matrix(m);
}

GrB_Matrix Graph_GetLabelMatrix(const Graph *g, int label_idx) {
	assert(g && label_idx < array_len(g->labels));
	RG_Matrix m = g->labels[label_idx];
	g->SynchronizeMatrix(g, m);
	return RG_Matrix_Get_GrB_Matrix(m);
}

GrB_Matrix Graph_GetRelationMatrix(const Graph *g, int relation_idx) {
	assert(g && (relation_idx == GRAPH_NO_RELATION || relation_idx < Graph_RelationTypeCount(g)));

	if(relation_idx == GRAPH_NO_RELATION) {
		return Graph_GetAdjacencyMatrix(g);
	} else {
		RG_Matrix m = g->relations[relation_idx];
		g->SynchronizeMatrix(g, m);
		return RG_Matrix_Get_GrB_Matrix(m);
	}
}

GrB_Matrix Graph_GetTransposedRelationMatrix(const Graph *g, int relation_idx) {
	assert(g && (relation_idx == GRAPH_NO_RELATION || relation_idx < Graph_RelationTypeCount(g)));

	if(relation_idx == GRAPH_NO_RELATION) {
		return Graph_GetTransposedAdjacencyMatrix(g);
	} else {
		assert(g->t_relations && "tried to retrieve nonexistent transposed matrix.");

		RG_Matrix m = g->t_relations[relation_idx];
		g->SynchronizeMatrix(g, m);
		return RG_Matrix_Get_GrB_Matrix(m);
	}
}

GrB_Matrix Graph_GetZeroMatrix(const Graph *g) {
	GrB_Index nvals;
	RG_Matrix z = g->_zero_matrix;
	g->SynchronizeMatrix(g, z);

	// Make sure zero matrix is indeed empty.
	GrB_Matrix grb_z = RG_Matrix_Get_GrB_Matrix(z);
	GrB_Matrix_nvals(&nvals, grb_z);
	assert(nvals == 0);
	return grb_z;
}

void Graph_Free(Graph *g) {
	assert(g);
	// Free matrices.
	Entity *en;
	DataBlockIterator *it;
	RG_Matrix_Free(g->_zero_matrix);
	RG_Matrix_Free(g->adjacency_matrix);
	RG_Matrix_Free(g->_t_adjacency_matrix);

    //printf("  Graph_Free: after Matrix Free\n"); fflush(stdout);

	_Graph_FreeRelationMatrices(g);
    //printf("  Graph_Free: after Free relation matrices\n"); fflush(stdout);
	array_free(g->relations);
	array_free(g->t_relations);
    //printf("  Graph_Free: after array free\n"); fflush(stdout);

	uint32_t labelCount = array_len(g->labels);
	for(int i = 0; i < labelCount; i++) {
		RG_Matrix_Free(g->labels[i]);
	}
	array_free(g->labels);
    //printf("  Graph_Free: after Free labels\n"); fflush(stdout);

	it = Graph_ScanNodes(g);
	while((en = (Entity *)DataBlockIterator_Next(it, NULL)) != NULL)
		FreeEntity(en);
    //printf("  Graph_Free: after Free entity\n"); fflush(stdout);

	DataBlockIterator_Free(it);
    //printf("  Graph_Free: after datablockiterator free\n"); fflush(stdout);

	it = Graph_ScanEdges(g);
	while((en = DataBlockIterator_Next(it, NULL)) != NULL)
		FreeEntity(en);
    //printf("  Graph_Free: after datablockiterator next\n"); fflush(stdout);

	DataBlockIterator_Free(it);
    //printf("  Graph_Free: after datablockiterator free #2\n"); fflush(stdout);

	// Free blocks.
	DataBlock_Free(g->nodes);
	DataBlock_Free(g->edges);
    //printf("  Graph_Free: after free blocks\n"); fflush(stdout);

	assert(pthread_mutex_destroy(&g->_writers_mutex) == 0);

	if(g->_writelocked) Graph_ReleaseLock(g);
	assert(pthread_rwlock_destroy(&g->_rwlock) == 0);
    //printf("  Graph_Free: after rwlock destroy\n"); fflush(stdout);

	rm_free(g);
    //printf("  Graph_Free: after rm free\n"); fflush(stdout);
}

