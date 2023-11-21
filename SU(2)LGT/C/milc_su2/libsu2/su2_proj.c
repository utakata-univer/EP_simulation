/*****************  su2_proj.c  (in su2.a) ******************************
*									*
* void su2_projector(su2_vector *a, su2_vector *b, su3_matrix *c);	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void su2_projector(su2_vector *a, su2_vector *b, su3_matrix *c) {
  register int i,j;

  for(i=0;i<2;i++)for(j=0;j<2;j++) { CMUL_J(a->c[i], b->c[j], c->e[i]) }
}
