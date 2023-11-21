/******************  s_m_s_vec.c  (in su2.a) ******************************
*									  *
* void scalar_mult_sub_su2_vector(su2_vector *a, su2_vector *b,           *
                                                 float s, su2_vector *c); *
* C <- A - s*B,  A, B and C vectors 						  *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_sub_su2_vector(su2_vector *a, su2_vector *b,
                                float s, su2_vector *c) {
   register int i;
   for(i=0;i<2;i++){
      c->c[i].real = a->c[i].real - s*b->c[i].real;
      c->c[i].imag = a->c[i].imag - s*b->c[i].imag;
   }
}
