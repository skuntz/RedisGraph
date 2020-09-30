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

     assert(NULL != iter);

     current_col = iter->current_col;
     const GrB_Index end_col = iter->end_col;
     if (end_col == current_col) {
          if (iter->k != iter->col_pattern_storage_len) {
               assert(iter->k == iter->col_pattern_len);
               return GrB_INDEX_OUT_OF_BOUNDS;
          }
          return GrB_SUCCESS;
     }

     const GrB_Matrix A = iter->A;
     const GrB_Index nr = iter->nr;
     const GrB_Index nc = iter->nc;

     info = GrB_Vector_new (&col, GrB_UINT64, nc);//nr);
     assert (GrB_SUCCESS == info);

     // Find the next non-empty column.
     for (; current_col < end_col; ++current_col) {
          info = GrB_extract (col, GrB_NULL, GrB_NULL, A, GrB_ALL, nr, current_col, GrB_DESC_T0);//GrB_NULL);
          if (info != GrB_SUCCESS) fprintf (stderr, "augh %d: %s\n", (int)info, GrB_error());
          assert (info == GrB_SUCCESS);
          info = GrB_Vector_nvals (&col_nv, col);
          assert (info == GrB_SUCCESS);
          if (col_nv) break;
     }

     iter->current_col = current_col;

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
          assert(NULL != new_col_pattern_storage);
          iter->col_pattern_storage = new_col_pattern_storage;
          iter->col_pattern_storage_len = col_nv;
     }
     iter->col_pattern = iter->col_pattern_storage;
     iter->col_pattern_len = col_nv;

     // Finally.  Extract the pattern.
     GrB_Index returned_nv = col_nv;
     info = GrB_Vector_extractTuples (iter->col_pattern, (bool*)GrB_NULL, &returned_nv, col);
     assert (GrB_SUCCESS == info);
     assert (returned_nv == col_nv);

     iter->k = 0;
     GrB_free (&col);
     return GrB_SUCCESS;
}

GrB_Info
iter_step (struct GBr_MatrixTupleIter_opaque* iter)
{
     assert (NULL != iter);

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
     assert (GrB_SUCCESS == info);
     info = GrB_Matrix_nrows (&iter->nr, iter->A);
     assert (iter->nr > 0);
     assert (GrB_SUCCESS == info);
     assert (begin_col <= iter->nc);
     assert (end_col <= iter->nc || end_col == IDX_ALL);
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
     assert (NULL != iter);
     GrB_Info info;
     memset (iter, 0, sizeof (*iter));
     info = iter_reuse (iter, A, begin_col, end_col);
     return info;
}

GrB_Info
iter_reset_range (struct GBr_MatrixTupleIter_opaque* iter, GrB_Index begin_col, GrB_Index end_col)
{
     assert (NULL != iter);
     assert (begin_col < iter->nc);
     assert (end_col <= iter->nc);
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
     assert (NULL != out);
     GrB_Info info = iter_init (out, A, 0, IDX_ALL);
     assert (GrB_SUCCESS == info);
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
     assert (NULL != iter);
     GrB_Info info;

     endColIdx += 1; // oy.
     return iter_reset_range (iter, startColIdx, endColIdx);
}

GrB_Info
GxB_MatrixTupleIter_iterate_range (GxB_MatrixTupleIter iter, 
                                   GrB_Index startColIdx, GrB_Index endColIdx)
{
     assert (NULL != iter);
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
     assert (NULL != iter);
     struct GBr_MatrixTupleIter_opaque *it = iter;
     if (iter_finishedp (iter)) { if (depleted) *depleted = true; return GrB_SUCCESS; }

     if (col) *col = it->current_col;
     assert(it->k < it->col_pattern_len);
     if (row) *row = it->col_pattern[it->k];
     return iter_step (iter);
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
     assert (NULL != iter);
     return iter_reuse (iter, A, 0, IDX_ALL);
}

