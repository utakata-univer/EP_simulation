/********************  s_m_hwvec.c  (in su2.a) ********************
*
* void scalar_mult_su2_hwvec(su2_half_wilson_vector *src, float s,
*	su2_half_wilson_vector *dest)
*  Multiply a half Wilson vector by a scalar
*  dest  <-  s*src
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_su2_hwvec(su2_half_wilson_vector *src, float s,
	su2_half_wilson_vector *dest) {
  register int i;
  for(i=0;i<2;i++) scalar_mult_su2_vector( &(src->h[i]), s, &(dest->h[i]));
}
