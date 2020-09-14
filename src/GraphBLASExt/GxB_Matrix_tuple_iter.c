//------------------------------------------------------------------------------
// GxB_MatrixTupleIter: Iterates over matrix none zero values
//------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdint.h>

// Should check that GrB_Index really is uint64_t
#define IDX_ALL ((GrB_Index)UINT64_MAX)

#include <assert.h>

#include "GraphBLAS.h"

#define SUBSYSTEM "iter: "

// XXX: Apparently _Generic and __VA_ARGS__ do not play well together.
#define ASSERT_ON_FAILURE(cond, str, errcode, ...) ASSERT_ON_FAILURE_(cond, __FILE__ ":" SUBSYSTEM str, errcode, do { __VA_ARGS__; } while (0))
#define ASSERT_ON_FAILURE_(cond, str, errcode, ...) do {bool augh = (cond); if (!augh) assert(NULL == str);  do { __VA_ARGS__; } while (0); return errcode; } while (0)

struct iter {
     GrB_Matrix A;
     GrB_Index nr, nc;
     GrB_Index begin_row, end_row;
     GrB_Index current_row, k;
     GrB_Index *row_pattern, row_pattern_len;
     GrB_Index *row_pattern_storage, row_pattern_storage_len;
};

typedef struct iter* GxB_MatrixTupleIter;

#define GxBr_NOT_OPAQUE
#include "GxB_Matrix_tuple_iter.h"

bool
iter_finishedp (const struct iter* iter)
{
     return iter->current_row == iter->end_row;
}

GrB_Info
iter_forward_row (struct iter* iter)
{
     GrB_Info info;
     GrB_Index current_row, row_nv = 0;
     GrB_Vector row;

     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);

     current_row = iter->current_row;
     const GrB_Index end_row = iter->end_row;
     if (end_row == current_row) {
          ASSERT_ON_FAILURE(iter->k == iter->row_pattern_len, "forwarding on non-empty row", 
                            GrB_INDEX_OUT_OF_BOUNDS);
          return GrB_SUCCESS;
     }

     const GrB_Matrix A = iter->A;
     const GrB_Index nc = iter->nc;

     info = GrB_Vector_new (&row, GrB_UINT64, iter->nc);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "new vector failed", info);

     // Find the next non-empty row.
     if (NULL == iter->row_pattern) { // Just beginning.
          assert (current_row == iter->begin_row);
          goto load_row;
     }
     do {
          ++current_row;
     load_row:
          info = GrB_extract (row, GrB_NULL, GrB_NULL, A, GrB_ALL, nc, current_row, GrB_DESC_T0);
          ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector", info, GrB_free (&row));
          info = GrB_Vector_nvals (&row_nv, row);
          ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector nvals", info, GrB_free (&row));
     } while (current_row < end_row && row_nv == 0);

     if (current_row == end_row) {
          // No more entries.
          GrB_free (&row);
          return GrB_SUCCESS;
     }

     // Ensure there's space.
     if (row_nv > iter->row_pattern_storage_len) {
          GrB_Index * new_row_pattern_storage;
          new_row_pattern_storage = realloc (iter->row_pattern_storage, 
                                             row_nv * sizeof (*iter->row_pattern_storage));
          ASSERT_ON_FAILURE(NULL != new_row_pattern_storage, 
                            "could not realloc row pattern storage", GrB_OUT_OF_MEMORY, GrB_free (&row));
          iter->row_pattern_storage = new_row_pattern_storage;
          iter->row_pattern_storage_len = row_nv;
     }
     iter->row_pattern = iter->row_pattern_storage;
     iter->row_pattern_len = row_nv;

     // Finally.  Extract the pattern.
     GrB_Index returned_nv = row_nv;
     info = GrB_Vector_extractTuples (iter->row_pattern, GrB_NULL, &returned_nv, row);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to extract vector pattern", info, GrB_free (&row));
     ASSERT_ON_FAILURE(returned_nv == row_nv, "vector nvals changed during execution", GrB_PANIC, GrB_free (&row));

     iter->k = 0;
     GrB_free (&row);
     return GrB_SUCCESS;
}

GrB_Info
iter_step (struct iter* iter)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);

     // Move forward within the row
     ++iter->k;
     assert (iter->k <= iter->row_pattern_len);
     
     if (iter->k == iter->row_pattern_len)
          return iter_forward_row (iter);

     return GrB_SUCCESS;
}

void
iter_val (GrB_Index* j, const struct iter* iter)
{
     // Really should check that the iterator is not null, finished, etc.
     // But apparently RedisGraph builds their release without NDEBUG, so performance here matters.
     *j = iter->row_pattern[iter->k];
}

GrB_Info
iter_reuse (struct iter* iter, GrB_Matrix A, GrB_Index begin_row, GrB_Index end_row)
{
     iter->row_pattern = NULL;
     iter->row_pattern_len = 0;
     iter->A = A;
     GrB_Info info = GrB_Matrix_nrows (&iter->nr, iter->A);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "nrows failed", info);
     info = GrB_Matrix_ncols (&iter->nc, iter->A);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "ncols failed", info);
     ASSERT_ON_FAILURE(begin_row < iter->nr, "starting iterator past the end", GrB_INVALID_INDEX);
     ASSERT_ON_FAILURE(end_row <= iter->nr || end_row == IDX_ALL, "ending the iterator past the end", GrB_INVALID_INDEX);
     iter->begin_row = begin_row;
     if (end_row == IDX_ALL)
          iter->end_row = iter->nr;
     else
          iter->end_row = end_row;
     return iter_forward_row (iter);
}

GrB_Info
iter_init (struct iter* iter, GrB_Matrix A, GrB_Index begin_row, GrB_Index end_row)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);
     GrB_Info info;
     memset (iter, 0, sizeof (*iter));
     return iter_reuse (iter, A, begin_row, end_row);
}

GrB_Info
iter_reset_range (struct iter* iter, GrB_Index begin_row, GrB_Index end_row)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);
     ASSERT_ON_FAILURE(begin_row < iter->nr, "starting iterator past the end", GrB_INVALID_INDEX);
     ASSERT_ON_FAILURE(end_row <= iter->nr, "ending the iterator past the end", GrB_INVALID_INDEX);
     iter->begin_row = begin_row;
     iter->current_row = begin_row;
     iter->end_row = end_row;
     iter->k = 0;
     iter->row_pattern = NULL;
     iter->row_pattern_len = 0;
     return iter_forward_row (iter);
}

GrB_Info
iter_free (struct iter* iter)
{
     ASSERT_ON_FAILURE(NULL != iter, "null iterator input", GrB_NULL_POINTER);
     if (iter->row_pattern_storage) free (iter->row_pattern_storage);
     memset (iter, 0, sizeof (*iter));
     return GrB_SUCCESS;
}


// Create a new iterator
GrB_Info
GxB_MatrixTupleIter_new (GxB_MatrixTupleIter *gxb_iter, GrB_Matrix A)
{
     struct iter *out;
     out = malloc (sizeof (*out));
     ASSERT_ON_FAILURE(out != NULL, "failed to allocate iterator", GrB_OUT_OF_MEMORY);
     GrB_Info info = iter_init (out, A, 0, IDX_ALL);
     ASSERT_ON_FAILURE(info == GrB_SUCCESS, "failed to initialize iterator", GrB_INVALID_OBJECT, free(out); *gxb_iter=NULL);
     *gxb_iter = out;
     return GrB_SUCCESS;
}

GrB_Info
GxB_MatrixTupleIter_free (GxB_MatrixTupleIter *gxb_iter)
{
     struct iter *it = *gxb_iter;
     if (NULL == gxb_iter) return GrB_SUCCESS;
     if (it->row_pattern_storage)
          free (it->row_pattern_storage);
     free (it);
     *gxb_iter = NULL;
     return GrB_SUCCESS;
}

static GrB_Info 
MatrixTupleIter_reset_ (GxB_MatrixTupleIter iter,
                       GrB_Index startRowIdx, GrB_Index endRowIdx)
{
     GrB_Info info;
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);

     endRowIdx += 1; // oy.
     return iter_reset_range (iter, startRowIdx, endRowIdx);
}

GrB_Info
GxB_MatrixTupleIter_iterate_range (GxB_MatrixTupleIter iter, 
                                   GrB_Index startRowIdx, GrB_Index endRowIdx)
{
     ASSERT_ON_FAILURE(NULL != iter, "NULL iterator", GrB_NULL_POINTER);
     return MatrixTupleIter_reset_ (iter, startRowIdx, endRowIdx);
}

GrB_Info 
GxB_MatrixTupleIter_iterate_row (GxB_MatrixTupleIter iter, GrB_Index rowIdx)
{
     return GxB_MatrixTupleIter_iterate_range (iter, rowIdx, rowIdx);
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
     struct iter *it = iter;
     if (iter_finishedp (iter)) { if (depleted) *depleted = true; return GrB_SUCCESS; }

     *row = it->current_row;
     assert(it->k < it->row_pattern_len);
     *col = it->row_pattern[it->k];
     return iter_forward_row (iter);
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

