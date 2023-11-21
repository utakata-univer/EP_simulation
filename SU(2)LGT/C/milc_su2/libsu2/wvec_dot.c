/******************  wvec_dot.c  (in su2.a) ******************************
*									*
* complex su2_wvec_dot(su2_wilson_vector *a, su2_wilson_vector *b)		*
* return dot product of two wilson_vectors					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

complex su2_wvec_dot(su2_wilson_vector *a, su2_wilson_vector *b) {
  complex temp1,temp2;
  register int i;
  temp1.real = temp1.imag = 0.0;
  for(i=0;i<4;i++){
    CMULJ_(a->d[i].c[0],b->d[i].c[0],temp2); CSUM(temp1,temp2);
    CMULJ_(a->d[i].c[1],b->d[i].c[1],temp2); CSUM(temp1,temp2);
  }
  return(temp1);
}
