/********************  s_m_wvec.c  (in su2.a) ********************
*
*void scalar_mult_su2_wvec(su2_wilson_vector *src, float s, 
*                                        su2_wilson_vector *dest)
*  Multiply a Wilson vector by a scalar
*  dest  <-  s*src
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_su2_wvec(su2_wilson_vector *src, float s, 
                                        su2_wilson_vector *dest) {
  register int i;
  for(i=0;i<4;i++) scalar_mult_su2_vector( &(src->d[i]), s, &(dest->d[i]));
}
