/*****************  wvec_rdot.c  (in su2.a) ******************************
*									*
* float su2_wvec_rdot(su2_wilson_vector *a, su2_wilson_vector *b)	        *
* return real part of dot product of two wilson_vectors			*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float su2_wvec_rdot(su2_wilson_vector *a, su2_wilson_vector *b) {
  register float temp1,temp2;
  register int i;
  temp2=0.0;
  for(i=0;i<4;i++){
    temp1 = a->d[i].c[0].real * b->d[i].c[0].real; temp2 += temp1;
    temp1 = a->d[i].c[0].imag * b->d[i].c[0].imag; temp2 += temp1;
    temp1 = a->d[i].c[1].real * b->d[i].c[1].real; temp2 += temp1;
    temp1 = a->d[i].c[1].imag * b->d[i].c[1].imag; temp2 += temp1;
  }
  return(temp2);
}
