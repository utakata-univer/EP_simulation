/**************  m_mat_hwvec.c  (in su2.a) ******************************
*									*
* void mult_mat_su2_hwvec(su3_matrix *mat,			        *
*	su2_half_wilson_vector *src, su2_half_wilson_vector *dest)      *
*  multiply a Wilson half-vector by a matrix				*
*  dest  <-  mat*src							*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_mat_su2_hwvec(su2_matrix *mat,
	su2_half_wilson_vector *src, su2_half_wilson_vector *dest) {

    mult_su2_mat_vec(mat, &(src->h[0]), &(dest->h[0]) );
    mult_su2_mat_vec(mat, &(src->h[1]), &(dest->h[1]) );
}
