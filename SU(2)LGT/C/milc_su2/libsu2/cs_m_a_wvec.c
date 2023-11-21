/********************  	cs_m_a_wvec.c  (in su2.a) ********************
*
*void c_scalar_mult_add_su2_wvec(su2_wilson_vector *src1, 
*                                su2_wilson_vector *src2,
* 	                         complex *s, wilson_vector *dest)
*  Multiply a Wilson vector by a complex scalar and add to another vector
* dest  <-  src1 + s*src2
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void c_scalar_mult_add_su2_wvec(su2_wilson_vector *src1, 
				su2_wilson_vector *src2,
				complex *phase, su2_wilson_vector *dest) {
   register int i,j;
   register complex t;
   for(i=0;i<4;i++){
      for(j=0;j<2;j++){
	 CMUL(src2->d[i].c[j], *phase, t);
	 CADD(src1->d[i].c[j], t, dest->d[i].c[j]);
      }
   }
}
