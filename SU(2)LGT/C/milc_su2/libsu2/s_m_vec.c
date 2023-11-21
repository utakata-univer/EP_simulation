/******************  s_m_vec.c  (in su2.a) ******************************
*									*
* void scalar_mult_su2_vector(su2_vector *a, float s, su2_vector *c);   *
* C <- s*A,  A and C vectors 						*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_su2_vector(su2_vector *a, float s, su2_vector *c) {
   register int i;
   for(i=0;i<2;i++){
	c->c[i].real = s*a->c[i].real;
	c->c[i].imag = s*a->c[i].imag;
    }
}
