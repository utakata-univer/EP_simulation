/*****************  m_amv_4dir.c  (in su2.a) *****************************
*									*
*  void mult_adj_su2_mat_vec_4dir(su2_matrix *mat, su2_vector *src,     *
*                                         su2_vector *dest);            *
*  Multiply an su3_vector by an array of four adjoint su2_matrices,	*
*  result in an array of four su2_vectors.				*
*  dest[i]  <-  A_adjoint[i] * src					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_adj_su2_mat_vec_4dir(su2_matrix *mat, su2_vector *src,
                                                     su2_vector *dest) {
    mult_adj_su2_mat_vec( mat+0, src, dest+0 );
    mult_adj_su2_mat_vec( mat+1, src, dest+1 );
    mult_adj_su2_mat_vec( mat+2, src, dest+2 );
    mult_adj_su2_mat_vec( mat+3, src, dest+3 );
}
