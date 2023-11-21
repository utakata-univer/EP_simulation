/******************  cs_m_a_vec.c  (in su2.a) ***************************
*									*
*  c_scalar_mult_add_su2vec():						*
*  multiply an su2 vector by a complex scalar and add it to another	*
*  vector:   v1 <- v1 + number*v2 					*
*  NOTE: non-standard order of target and source v1 and v2!             *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void c_scalar_mult_add_su2vec(su2_vector *v1, complex *phase, su2_vector *v2) {
  register int i;
  complex t;

  for(i=0;i<2;i++){
    CMUL(v2->c[i], *phase, t);
    CADD(v1->c[i], t, v1->c[i]);
  }
}
