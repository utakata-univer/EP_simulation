/*****************  su3mat_copy.c  (in su3.a) ***************************
*									*
* void su3mat_copy(a,b) su3_matrix *a,*b;				*
* Copy an su3 matrix:  B <- A   						*
*/
#include "complex.h"
#include "su3.h"

/* Copy a su3 matrix:  b <- a   */
void su3mat_copy(a,b) su3_matrix *a,*b; {
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	b->e[i][j].real = a->e[i][j].real;
	b->e[i][j].imag = a->e[i][j].imag;
    }
}
