/******************  s_m_sum_vec.c  (in su2.a) ******************************
*									    *
* void scalar_mult_sum_su2_vector(su2_vector *a, su2_vector *b, float s);   *
* A <- A + s*B,  A and C vectors 					    *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_sum_su2_vector(su2_vector *a, su2_vector *b, float s) {
   register int i;
   for(i=0;i<2;i++){
      a->c[i].real += s*b->c[i].real;
      a->c[i].imag += s*b->c[i].imag;
   }
}
