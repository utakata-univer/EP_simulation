/*****************  s_m_atm_wvec.c  (in su2.a) ********************
*
*void scalar_mult_addtm_su2_wvec(su2_wilson_vector *src1, 
*                                su2_wilson_vector *src2,
*                                float s, su2_wilson_vector *dest)
*  Multiply a Wilson vector by a scalar and add to minus one times
*   another vector
* dest  <-  (-1)*src1 + s*src2
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_addtm_su2_wvec(su2_wilson_vector *src1, 
				su2_wilson_vector *src2,
				float s, su2_wilson_vector *dest) {
  register int i,j;
  for(i=0;i<4;i++){	/*spins*/
    for(j=0;j<2;j++){  /*colors*/
      dest->d[i].c[j].real = -src1->d[i].c[j].real +
	                        s*src2->d[i].c[j].real;
      dest->d[i].c[j].imag = -src1->d[i].c[j].imag +
	                        s*src2->d[i].c[j].imag;
    }
  }
}
