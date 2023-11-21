/********************  subvec.c  (in su2.a) *****************************
*									*
*  Subtract two SU2 vectors							*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void sub_su2_vector(su2_vector *a, su2_vector *b, su2_vector *c) {
  register int i;
  for(i=0;i<2;i++){ CSUB(a->c[i], b->c[i], c->c[i]); }
}
