/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/

#include "bulk_insert.h"
#include "RG.h"
#include "../schema/schema.h"
#include "../util/rmalloc.h"
#include "../datatypes/array.h"
#include <errno.h>

#include <sys/time.h>
int timeval_subtract(struct timeval *result,
struct timeval end,
struct timeval start)
{
if (start.tv_usec < end.tv_usec) {
int nsec = (end.tv_usec - start.tv_usec) / 1000000 + 1;
end.tv_usec -= 1000000 * nsec;
end.tv_sec += nsec;
}
if (start.tv_usec - end.tv_usec > 1000000) {
int nsec = (end.tv_usec - start.tv_usec) / 1000000;
end.tv_usec += 1000000 * nsec;
end.tv_sec -= nsec;
}

result->tv_sec = end.tv_sec - start.tv_sec;
result->tv_usec = end.tv_usec - start.tv_usec;

return end.tv_sec < start.tv_sec;
}

void set_exec_time(int end)
{
    static struct timeval time_start;
    struct timeval time_end;
    struct timeval time_diff;

    if (end) {
        gettimeofday(&time_end, NULL);
        if (timeval_subtract(&time_diff, time_end, time_start) == 0) {
            if (end == 1)
                printf("\nexec time: %1.2fs\n",
                       time_diff.tv_sec + (time_diff.tv_usec / 1000000.0f));
            else if (end == 2)
                printf("%1.2fs",
                       time_diff.tv_sec + (time_diff.tv_usec / 1000000.0f));
        }
        return;
    }
    gettimeofday(&time_start, NULL);
}

void start_exec_timer()
{
    set_exec_time(0);
}

void print_exec_timer()
{
    set_exec_time(1);
}





// The first byte of each property in the binary stream
// is used to indicate the type of the subsequent SIValue
typedef enum {
	BI_NULL = 0,
	BI_BOOL = 1,
	BI_DOUBLE = 2,
	BI_STRING = 3,
	BI_LONG = 4,
	BI_ARRAY = 5,
} TYPE;

// Read the header of a data stream to parse its property keys and update schemas.
static Attribute_ID *_BulkInsert_ReadHeader(GraphContext *gc, SchemaType t,
											const char *data, size_t *data_idx,
											int *label_id, unsigned int *prop_count) {
	/* Binary header format:
	 * - entity name : null-terminated C string
	 * - property count : 4-byte unsigned integer
	 * [0..property_count] : null-terminated C string
	 */

    int msec;
    printf("_BulkInsert_ReadHeader: \n"); fflush(stdout);
    start_exec_timer();

	// First sequence is entity name
	const char *name = data + *data_idx;
	*data_idx += strlen(name) + 1;
	Schema *schema = GraphContext_GetSchema(gc, name, t);
	if(schema == NULL) schema = GraphContext_AddSchema(gc, name, t);
	*label_id = schema->id;

	// Next 4 bytes are property count
	*prop_count = *(unsigned int *)&data[*data_idx];
	*data_idx += sizeof(unsigned int);

	if(*prop_count == 0) return NULL;
	Attribute_ID *prop_indicies = malloc(*prop_count * sizeof(Attribute_ID));

	// The rest of the line is [char *prop_key] * prop_count
	for(unsigned int j = 0; j < *prop_count; j ++) {
		char *prop_key = (char *)data + *data_idx;
		*data_idx += strlen(prop_key) + 1;

		// Add properties to schemas
		prop_indicies[j] = GraphContext_FindOrAddAttribute(gc, prop_key);
	}
    print_exec_timer();

	return prop_indicies;
}

// Read an SIValue from the data stream and update the index appropriately
static inline SIValue _BulkInsert_ReadProperty(const char *data, size_t *data_idx) {
	/* Binary property format:
	 * - property type : 1-byte integer corresponding to TYPE enum
	 * - Nothing if type is NULL
	 * - 1-byte true/false if type is boolean
	 * - 8-byte double if type is double
	 * - 8-byte integer if type is integer
	 * - Null-terminated C string if type is string
	 * - 8-byte array length followed by N values if type is array
	 */

    int msec;
    printf("_BulkInsert_ReadProperty: \n"); fflush(stdout);
    start_exec_timer();

	SIValue v = SI_NullVal();
	TYPE t = data[*data_idx];
	*data_idx += 1;

	if(t == BI_NULL) {
		v = SI_NullVal();
	} else if(t == BI_BOOL) {
		bool b = data[*data_idx];
		*data_idx += 1;
		v = SI_BoolVal(b);
	} else if(t == BI_DOUBLE) {
		double d = *(double *)&data[*data_idx];
		*data_idx += sizeof(double);
		v = SI_DoubleVal(d);
	} else if(t == BI_LONG) {
		int64_t d = *(int64_t *)&data[*data_idx];
		*data_idx += sizeof(int64_t);
		v = SI_LongVal(d);
	} else if(t == BI_STRING) {
		const char *s = data + *data_idx;
		*data_idx += strlen(s) + 1;
		// The string itself will be cloned when added to the GraphEntity properties.
		v = SI_ConstStringVal((char *)s);
	} else if(t == BI_ARRAY) {
		// The first 8 bytes of a received array will be the array length.
		int64_t len = *(int64_t *)&data[*data_idx];
		*data_idx += sizeof(int64_t);
		v = SIArray_New(len);
		for(uint i = 0; i < len; i ++) {
			// Convert every element and add to array.
			SIArray_Append(&v, _BulkInsert_ReadProperty(data, data_idx));
		}
	} else {
		ASSERT(false);
	}
	print_exec_timer();
	return v;
}

int _BulkInsert_ProcessNodeFile(RedisModuleCtx *ctx, GraphContext *gc, const char *data,
								size_t data_len) {

    int msec;
    printf("_BulkInsert_ProcessNodeFile: \n"); fflush(stdout);
    start_exec_timer();

	size_t data_idx = 0;
        Graph *g = gc->g;

	int label_id;
	unsigned int prop_count;
	Attribute_ID *prop_indicies = _BulkInsert_ReadHeader(gc, SCHEMA_NODE, data, &data_idx, &label_id,
														 &prop_count);
        const size_t post_header_data_idx = data_idx;

        /* Two passes.  The first counts the number of vertices in *this*
           file/chunk.  The second collects the properties. */
        uint64_t n_to_alloc = 0;
	while(data_idx < data_len) {
             for(unsigned int i = 0; i < prop_count; i++) _BulkInsert_ReadProperty(data, &data_idx);
             ++n_to_alloc;
	}

        const uint64_t start_id = Graph_BulkCreateNodes (g, label_id, n_to_alloc);

        /* Pre-allocate entity space.  Can shrink if NULLs appear. */
        if (prop_count > 0) {
             for (uint64_t k = start_id; k < start_id + n_to_alloc; ++k) {
                  Entity *en = DataBlock_GetItem (g->nodes, k);
                  ASSERT (en);
                  en->properties = rm_malloc (prop_count * sizeof (EntityProperty));
                  ASSERT (en->properties);
             }

             uint64_t id = start_id;
             data_idx = post_header_data_idx;
             while (data_idx < data_len) {
                  Entity *en = DataBlock_GetItem (g->nodes, id);
                  ASSERT (en);
                  ASSERT (en->prop_count == 0); // Allocated above, so should be zero.
                  int prop_idx = 0;
                  for(unsigned int i = 0; i < prop_count; i++) {
                       SIValue value = _BulkInsert_ReadProperty(data, &data_idx);
                       // Cypher does not support NULL as a property value.
                       // If we encounter one here, simply skip it.
                       if(SI_TYPE(value) == T_NULL) continue;
                       en->properties[prop_idx].id = prop_indicies[i];
                       en->properties[prop_idx].value = SI_CloneValue(value);
                       ++prop_idx;
                  }
                  en->prop_count = prop_idx;
                  if (prop_idx != prop_count) {
                       // Shrink-wrap the allocation.
                       EntityProperty *tmp = rm_realloc (en->properties, prop_idx * sizeof (*tmp));
                       ASSERT(tmp != NULL);
                       if (tmp != NULL)
                            en->properties = tmp; /* A "soft" failure.  The memory is there but wasted. */
                  }
                  ++id;
             }
        }

	free(prop_indicies);
    print_exec_timer();
	return BULK_OK;
}

int _BulkInsert_ProcessRelationFile(RedisModuleCtx *ctx, GraphContext *gc, const char *data,
                                    size_t data_len, NodeID* I, NodeID *J) {

    int msec;
    printf("_BulkInsert_ProcessRelationFile: \n"); fflush(stdout);
    start_exec_timer();

	size_t data_idx = 0;

	int reltype_id;
	unsigned int prop_count;
	// Read property keys from header and update schema
	Attribute_ID *prop_indicies = _BulkInsert_ReadHeader(gc, SCHEMA_EDGE, data, &data_idx, &reltype_id,
														 &prop_count);
        const size_t post_header_data_idx = data_idx;

        uint64_t n_relations_in_chunk = 0;

        if (prop_count > 0) {
             while (data_idx < data_len) {
                  ++n_relations_in_chunk;
                  data_idx += 2 * sizeof (NodeID); // Skip vertices.
                  for(unsigned int i = 0; i < prop_count; i ++)
			_BulkInsert_ReadProperty(data, &data_idx);
             }
        } else {
             n_relations_in_chunk = (data_len - post_header_data_idx) / (2 * sizeof (NodeID));
        }
        if (n_relations_in_chunk == 0) return BULK_OK;

        size_t k = 0;
        data_idx = post_header_data_idx;
        while (data_idx < data_len) {
             I[k] = *(NodeID *)&data[data_idx];
             data_idx += sizeof(NodeID);
             J[k] = *(NodeID *)&data[data_idx];
             data_idx += sizeof(NodeID);
             ++k;

             if (prop_count > 0)
                  for(unsigned int i = 0; i < prop_count; i ++)
			_BulkInsert_ReadProperty(data, &data_idx);
        }

        EdgeID start_id = Graph_BulkConnectNodes (gc->g, I, J, n_relations_in_chunk, reltype_id);
        // If -1, didn't change anything.  Probably should trigger garbage collection.
        if (start_id == (EdgeID)-1) return BULK_OK;

        /* Pre-allocate entity space.  Can shrink if NULLs appear. */
        if (prop_count > 0) {
             for (uint64_t id = start_id; id < start_id + n_relations_in_chunk; ++id) {
                  Entity *en = DataBlock_GetItem (gc->g->edges, id);
                  if (en) // Some may not have been "created" if not a multi-graph.
                       en->properties = rm_malloc (prop_count * sizeof (EntityProperty));
             }

             uint64_t id = start_id;
             data_idx = post_header_data_idx;
             while (data_idx < data_len) {
                  data_idx += 2 * sizeof (NodeID);
                  Entity *en = DataBlock_GetItem (gc->g->edges, id);
                  /* en could be NULL if this edge was deleted in a
                     non-multi-graph.  If non-NULL, the properties are
                     uninitialized, and prop_count == 0.  The
                     properties still need "read" from data[] because
                     a later edge may not be deleted.
                  */
                  ASSERT (en == NULL || en->prop_count == 0);
                  int prop_idx = 0;
                  for(unsigned int i = 0; i < prop_count; i++) {
                       SIValue value = _BulkInsert_ReadProperty(data, &data_idx);
                       // Cypher does not support NULL as a property value.
                       // If we encounter one here, simply skip it.
                       if(SI_TYPE(value) == T_NULL) continue;
                       if (en) {
                            en->properties[prop_idx].id = prop_indicies[i];
                            en->properties[prop_idx].value = SI_CloneValue(value);
                       }
                       ++prop_idx;
                  }
                  if (en) {
                       en->prop_count = prop_idx;
                       if (prop_idx != prop_count) {
                            if (en->properties) {
                                 // Shrink-wrap the allocation.
                                 EntityProperty *tmp = rm_realloc (en->properties, prop_idx * sizeof (*tmp));
                                 ASSERT(tmp != NULL);
                                 if (tmp != NULL)
                                      en->properties = tmp; /* A "soft" failure.  The memory is there but wasted. */
                            } else {
                                 ASSERT(prop_idx == 0);
                                 ASSERT(prop_count == 0);
                            }
                       }
                  }
                  ++id;
             }
        }

	free(prop_indicies);
    print_exec_timer();
	return BULK_OK;
}

int _BulkInsert_InsertNodes(RedisModuleCtx *ctx, GraphContext *gc, int token_count,
							RedisModuleString ***argv, int *argc) {
    int msec;
    printf("_BulkInsert_InsertNodes: \n"); fflush(stdout);
    start_exec_timer();

	int rc;
	for(int i = 0; i < token_count; i ++) {
		size_t len;
		// Retrieve a pointer to the next binary stream and record its length
		const char *data = RedisModule_StringPtrLen(**argv, &len);
		*argv += 1;
		*argc -= 1;
		rc = _BulkInsert_ProcessNodeFile(ctx, gc, data, len);
		UNUSED(rc);
		ASSERT(rc == BULK_OK);
	}
	print_exec_timer();
	return BULK_OK;
}

int _BulkInsert_Insert_Edges(RedisModuleCtx *ctx, GraphContext *gc, int token_count,
                             RedisModuleString ***argv, int *argc,
                             NodeID* I, NodeID* J) {

    int msec;
    printf("_BulkInsert_Insert_Edges: \n"); fflush(stdout);
    start_exec_timer();

	int rc;
	for(int i = 0; i < token_count; i ++) {
		size_t len;
		// Retrieve a pointer to the next binary stream and record its length
		const char *data = RedisModule_StringPtrLen(**argv, &len);
		*argv += 1;
		*argc -= 1;
		rc = _BulkInsert_ProcessRelationFile(ctx, gc, data, len, I, J);
		UNUSED(rc);
		ASSERT(rc == BULK_OK);
	}
	print_exec_timer();
	return BULK_OK;
}

int BulkInsert(RedisModuleCtx *ctx, GraphContext *gc, RedisModuleString **argv, int argc,
               long long nodes_in_query, long long relations_in_query) {

	if(argc < 2) {
		RedisModule_ReplyWithError(ctx, "Bulk insert format error, failed to parse bulk insert sections.");
		return BULK_FAIL;
	}

    time_t start, start1, end;
	int msec;
	start = clock();
	printf("BulkInsert: start = %ld \n", start); fflush(stdout);

	// Read the number of node tokens
	long long node_token_count;
	long long relation_token_count;
	if(RedisModule_StringToLongLong(*argv++, &node_token_count)  != REDISMODULE_OK) {
		RedisModule_ReplyWithError(ctx, "Error parsing number of node descriptor tokens.");
		return BULK_FAIL;
	}

	if(RedisModule_StringToLongLong(*argv++, &relation_token_count)  != REDISMODULE_OK) {
		RedisModule_ReplyWithError(ctx, "Error parsing number of relation descriptor tokens.");
		return BULK_FAIL;
	}
	argc -= 2;

	if(node_token_count > 0) {
	    start_exec_timer();
	    start1 = clock();
	    printf("BulkInsert_InsertNodes: \n"); fflush(stdout);
		int rc = _BulkInsert_InsertNodes(ctx, gc, node_token_count, &argv, &argc);
        end = clock();
        msec = (int) (end - start1) * 1000 / CLOCKS_PER_SEC;
        printf(" total time %d sec %d ms\n", msec/1000, msec % 1000); fflush(stdout);
        print_exec_timer();
		if(rc != BULK_OK) {
			return BULK_FAIL;
		} else if(argc == 0) {
			return BULK_OK;
		}

	}

	if(relation_token_count > 0) {
        start_exec_timer();
	    start1 = clock();
        printf("BulkInsert_InsertEdges: \n"); fflush(stdout);

                NodeID *I, *J;  // Imported vertex scratch space.
                I = rm_malloc (2 * relations_in_query * sizeof (*I));
                if (!I) return BULK_FAIL;
                J = &I[relations_in_query];
		int rc = _BulkInsert_Insert_Edges(ctx, gc, relation_token_count, &argv, &argc, I, J);
        end = clock();
        msec = (int) (end - start1) * 1000 / CLOCKS_PER_SEC;
        printf("   total time %d sec %d ms\n", msec/1000, msec % 1000); fflush(stdout);
        print_exec_timer();
		if(rc != BULK_OK) {
			return BULK_FAIL;
		} else if(argc == 0) {
			return BULK_OK;
		}


	}

	ASSERT(argc == 0);



	return BULK_OK;
}

