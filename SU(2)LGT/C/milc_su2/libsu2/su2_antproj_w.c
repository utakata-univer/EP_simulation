/*****************  su2_antproj_w.c  (in su2.a) ******************************
*									*
* void su2_antiherm_projector_w(su2_wilson_vector *a, 
*                             su2_wilson_vector *b, su2_matrix *c);	*
* C  <- sum over spins of outer product of A.d[i] and B.d[i]		*
*  C_ij = sum( A_i * B_adjoint_j )                                      *
* This version projects onto the traceless antihermitian part		*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void su2_antiherm_projector_w(su2_wilson_vector *a, 
                             su2_wilson_vector *b, su2_matrix *c) {
  register int k;
  su2_matrix cc;
  for(k=0;k<4;k++)c->e[k]=0.0;
  for(k=0;k<4;k++){
      su2_antiherm_projector(&(a->d[k]), &(b->d[k]), &cc ); 
      add_su2_matrix ( c, &cc, c );
  }
}

