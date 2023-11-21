/******************  magsq_su2vec.c  (in su2.a) *************************
*                                                                       *
* float magsq_su2vec(su2_vector *a);                                    *
* return squared magnitude of an SU2 vector
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float magsq_su2vec(su2_vector *a) {
   register float sum;
   register int i;
   for(i=0,sum=0.0;i<2;i++) 
      sum += a->c[i].real*a->c[i].real + a->c[i].imag*a->c[i].imag;
   return(sum);
}

