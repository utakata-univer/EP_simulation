/******************  su2_adjoint.c  (in su2.a) **************************
*                                                                       *
* void su2_adjoint(a,b) su2_matrix *a,*b;                               *
* B  <- A_adjoint,  adjoint of an SU2 matrix                            *
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

/* adjoint of an SU2 matrix */
void su2_adjoint(su2_matrix *a, su2_matrix *b) {
   register int i;

   b->e[0] = a->e[0]; 
   for(i=1;i<=3;i++) b->e[i] = -a->e[i];
}
