/********************  s_m_a_wvec.c  (in su2.a) ********************
*
*void scalar_mult_add_su2_wvec(su2_wilson_vector *src1, 
*       su2_wilson_vector *src2, float s, su2_wilson_vector *dest)
*  Multiply a Wilson vector by a scalar and add to another vector
*  dest  <-  src1 + s*src2
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"


void scalar_mult_add_su2_wvec(su2_wilson_vector *src1, 
			      su2_wilson_vector *src2,
			      float s, su2_wilson_vector *dest) {
  register int i;
  for(i=0;i<4;i++) scalar_mult_add_su2_vector( &(src1->d[i]), &(src2->d[i]),
	s, &(dest->d[i]));
}
