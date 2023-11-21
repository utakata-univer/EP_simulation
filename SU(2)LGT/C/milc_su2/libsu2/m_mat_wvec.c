/******************  m_mat_wvec.c  (in su2.a) ********************
*									*
*void mult_mat_su2_wilson_vec(su2_matrix *mat, su2_wilson_vector *src,  *
                                           su2_wilson_vector *dest)	*
*  multiply a Wilson vector by a matrix					*
* dest  <-  mat*src							*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_mat_su2_wilson_vec(su2_matrix *mat, su2_wilson_vector *src,
                                          su2_wilson_vector *dest) {
  register int i;
    for(i=0;i<4;i++) mult_su2_mat_vec(mat, &(src->d[i]), &(dest->d[i]) );
}
