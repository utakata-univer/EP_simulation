/******************  realtr.c  (in su2.a) *******************************
*									*
* float realtrace_su2(su2_matrix *a, su2_matrix *b)			*
* return Re( Tr( A_adjoint * B )  					*
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float realtrace_su2(su2_matrix *a, su2_matrix *b) {
  return( 2.*(a->e[0]*b->e[0] + 
            (a->e[1]*b->e[1] + a->e[2]*b->e[2] + a->e[3]*b->e[3])));
}


