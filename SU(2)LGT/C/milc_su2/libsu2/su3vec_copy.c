/*****************  su3vec_copy.c  (in su3.a) ***************************
*									*
* void su3vec_copy(a,b) su3_vector *a,*b;				*
* Copy an su3 vector:  B <- A   					*
*/
#include "complex.h"
#include "su3.h"

/* Copy a su3 vector:  b <- a   */
void su3vec_copy(a,b) su3_vector *a,*b; {
register int i;
    for(i=0;i<3;i++){
	b->c[i].real = a->c[i].real;
	b->c[i].imag = a->c[i].imag;
    }
}
