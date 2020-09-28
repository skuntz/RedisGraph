//------------------------------------------------------------------------------
// GxB_MatrixTupleIter: Iterates over matrix none zero values
//------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdint.h>

#include <stdio.h>

#include <assert.h>

#include "GraphBLAS.h"

// Not great, but handy.
#define IDX_ALL ((uint64_t)-1)

#define SUBSYSTEM "iter: "

// XXX: Apparently _Generic and __VA_ARGS__ do not play well together.
#if 1
#define ASSERT_ON_FAILURE(cond, str, errcode, ...) ASSERT_ON_FAILURE_(cond, __FILE__ ":" SUBSYSTEM str, errcode, do { __VA_ARGS__; } while (0))
#define ASSERT_ON_FAILURE_(cond, str, errcode, ...) do {bool augh = (cond); if (!augh) assert(NULL == str);  do { __VA_ARGS__; } while (0); return errcode; } while (0)
#else
#define ASSERT_ON_FAILURE(...)
#endif

struct GBr_MatrixTupleIter_opaque {
     GrB_Matrix A;
     GrB_Index nr, nc;
     GrB_Index begin_col, end_col;
     GrB_Index current_col, k;
     GrB_Index *col_pattern, col_pattern_len;
     GrB_Index *col_pattern_storage, col_pattern_storage_len;
};

#include "GxB_Matrix_tuple_iter.h"

bool
iter_finishedp (const struct GBr_MatrixTupleIter_opaque* iter)
{
     return iter->current_col == iter->end_col;
}

GrB_Info
iter_forward_col (struct GBr_MatrixTupleIter_opaque* iter)
{
     GrB_Info info;
     GrB_Index current_col, col_nv = 0;
     GrB_Vector col;

     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);

     current_col = iter->current_col;
     const GrB_Index end_col = iter->end_col;
     if (end_col == current_col) {
          ASSERT_ON_FAILURE(iter->k == iter->col_pattern_len, "forwarding on non-empty column", 
                            GrB_INDEX_OUT_OF_BOUNDS);
          return GrB_SUCCESS;
     }

     const GrB_Matrix A = iter->A;
     const GrB_Index nr = iter->nr;
     const GrB_Index nc = iter->nc;

     info = GrB_Vector_new (&col, GrB_UINT64, nr);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "new vector failed", info);

     // Find the next non-empty column.
     if (NULL == iter->col_pattern) { // Just beginning.
          assert (current_col == iter->begin_col);
          goto load_col;
     }
     do {
          ++current_col;
     load_col:
          info = GrB_extract (col, GrB_NULL, GrB_NULL, A, GrB_ALL, nr, current_col, GrB_NULL);
          ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector", info, GrB_free (&col));
          info = GrB_Vector_nvals (&col_nv, col);
          ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector nvals", info, GrB_free (&col));
     } while (current_col < end_col && col_nv == 0);

     if (current_col == end_col) {
          // No more entries.
          GrB_free (&col);
          return GrB_SUCCESS;
     }

     // Ensure there's space.
     if (col_nv > iter->col_pattern_storage_len) {
          GrB_Index * new_col_pattern_storage;
          new_col_pattern_storage = realloc (iter->col_pattern_storage, 
                                             col_nv * sizeof (*iter->col_pattern_storage));
          ASSERT_ON_FAILURE(NULL != new_col_pattern_storage, 
                            "could not realloc col pattern storage", GrB_OUT_OF_MEMORY, GrB_free (&col));
          iter->col_pattern_storage = new_col_pattern_storage;
          iter->col_pattern_storage_len = col_nv;
     }
     iter->col_pattern = iter->col_pattern_storage;
     iter->col_pattern_len = col_nv;

     // Finally.  Extract the pattern.
     GrB_Index returned_nv = col_nv;
     info = GrB_Vector_extractTuples (GrB_NULL, iter->col_pattern, &returned_nv, col);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector pattern", info, GrB_free (&col));
     ASSERT_ON_FAILURE(returned_nv == col_nv, "vector nvals changed during execution", GrB_PANIC, GrB_free (&col));

     iter->k = 0;
     GrB_free (&col);
     return GrB_SUCCESS;
}

GrB_Info
iter_step (struct GBr_MatrixTupleIter_opaque* iter)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);

     // Move forward within the column
     ++iter->k;
     assert (iter->k <= iter->col_pattern_len);
     
     if (iter->k == iter->col_pattern_len)
          return iter_forward_col (iter);

     return GrB_SUCCESS;
}

void
iter_val (GrB_Index* j, const struct GBr_MatrixTupleIter_opaque* iter)
{
     // Really should check that the iterator is not null, finished, etc.
     // But apparently RedisGraph builds their release without NDEBUG, so performance here matters.
     *j = iter->col_pattern[iter->k];
}

GrB_Info
iter_reuse (struct GBr_MatrixTupleIter_opaque* iter, GrB_Matrix A, GrB_Index begin_col, GrB_Index end_col)
{
     iter->col_pattern = NULL;
     iter->col_pattern_len = 0;
     iter->A = A;
     GrB_Info info = GrB_Matrix_ncols (&iter->nc, iter->A);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "ncols failed", info);
     info = GrB_Matrix_nrows (&iter->nr, iter->A);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "nrows failed", info);
     ASSERT_ON_FAILURE(begin_col <= iter->nc, "starting iterator past the end", GrB_INVALID_INDEX);
     ASSERT_ON_FAILURE(end_col <= iter->nc || end_col == IDX_ALL, "ending the iterator past the end", GrB_INVALID_INDEX);
     iter->begin_col = begin_col;
     if (end_col == IDX_ALL)
          iter->end_col = iter->nc;
     else
          iter->end_col = end_col;
     if (iter->nc)
          return iter_forward_col (iter);
     else
          return GrB_SUCCESS;
}

GrB_Info
iter_init (struct GBr_MatrixTupleIter_opaque* iter, GrB_Matrix A, GrB_Index begin_col, GrB_Index end_col)
{
     //ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);
     GrB_Info info;
     memset (iter, 0, sizeof (*iter));
     info = iter_reuse (iter, A, begin_col, end_col);
     return info;
}

GrB_Info
iter_reset_range (struct GBr_MatrixTupleIter_opaque* iter, GrB_Index begin_col, GrB_Index end_col)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);
     ASSERT_ON_FAILURE(begin_col < iter->nc, "starting iterator past the end", GrB_INVALID_INDEX);
     ASSERT_ON_FAILURE(end_col <= iter->nc, "ending the iterator past the end", GrB_INVALID_INDEX);
     iter->begin_col = begin_col;
     iter->current_col = begin_col;
     iter->end_col = end_col;
     iter->k = 0;
     iter->col_pattern = NULL;
     iter->col_pattern_len = 0;
     return iter_forward_col (iter);
}

// Create a new iterator
GrB_Info
GxB_MatrixTupleIter_new (GxB_MatrixTupleIter *gxb_iter, GrB_Matrix A)
{
     if (gxb_iter == NULL) { return GrB_NULL_POINTER; }
     struct GBr_MatrixTupleIter_opaque *out;
     out = malloc (sizeof (*out));
     //ASSERT_ON_FAILURE(out != NULL, "failed to allocate iterator", GrB_OUT_OF_MEMORY);
     GrB_Info info = iter_init (out, A, 0, IDX_ALL);
     //ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to initialize iterator", GrB_INVALID_OBJECT, free(out); *gxb_iter=NULL);
     assert (out != NULL);
     *gxb_iter = out;
     assert (*gxb_iter != NULL);
     return GrB_SUCCESS;
}

GrB_Info
GxB_MatrixTupleIter_free (GxB_MatrixTupleIter *gxb_iter)
{
     struct GBr_MatrixTupleIter_opaque *it = *gxb_iter;
     if (NULL == gxb_iter) return GrB_SUCCESS;
     if (it->col_pattern_storage)
          free (it->col_pattern_storage);
     free (it);
     *gxb_iter = NULL;
     return GrB_SUCCESS;
}

GrB_Info 
MatrixTupleIter_reset_ (GxB_MatrixTupleIter iter,
                        GrB_Index startColIdx, GrB_Index endColIdx)
{
     GrB_Info info;
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);

     endColIdx += 1; // oy.
     return iter_reset_range (iter, startColIdx, endColIdx);
}

GrB_Info
GxB_MatrixTupleIter_iterate_range (GxB_MatrixTupleIter iter, 
                                   GrB_Index startColIdx, GrB_Index endColIdx)
{
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);
     return MatrixTupleIter_reset_ (iter, startColIdx, endColIdx);
}

// Advance iterator
GrB_Info GxB_MatrixTupleIter_next
(
	GxB_MatrixTupleIter iter,      // iterator to consume
	GrB_Index *row,                 // optional output row index
	GrB_Index *col,                 // optional output column index
	bool *depleted                  // indicate if iterator depleted
) {
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);
     struct GBr_MatrixTupleIter_opaque *it = iter;
     if (iter_finishedp (iter)) { if (depleted) *depleted = true; return GrB_SUCCESS; }

     if (col) *col = it->current_col;
     assert(it->k < it->col_pattern_len);
     if (row) *row = it->col_pattern[it->k];
     return iter_forward_col (iter);
}

// Reset iterator
GrB_Info GxB_MatrixTupleIter_reset
(
	GxB_MatrixTupleIter iter       // iterator to reset
) {
     return MatrixTupleIter_reset_ (iter, 0, IDX_ALL);
}

// Update iterator to scan given matrix
GrB_Info GxB_MatrixTupleIter_reuse
(
	GxB_MatrixTupleIter iter,      // iterator to update
	GrB_Matrix A                    // matrix to scan
) {
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);
     return iter_reuse (iter, A, 0, IDX_ALL);
}

