/*****************  su2_rdot.c  (in su2.a) ******************************
*                                                                       *
* float su2_rdot(su2_vector *a, su2_vector *b);                         *
* return real part of dot product of two su2_vectors                    *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float su2_rdot(su2_vector *a, su2_vector *b) {
   register float sum;
   register int i;
   for(i=0,sum=0.0; i<2; i++) 
      sum += a->c[i].real*b->c[i].real + a->c[i].imag*b->c[i].imag;
   return(sum);
}
