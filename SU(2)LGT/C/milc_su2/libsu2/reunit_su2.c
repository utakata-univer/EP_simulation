/******************  reunit_su2.c  (in su2.a) ***************************
*									*
* void reunit_su2( su2_matrix *a )					*
* reunitarize su2_matrix so that det a = 1
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"
double sqrt(double);

void reunit_su2(su2_matrix *a) {
  register int i;
  float asq=0.;

  for(i=0; i<4; ++i) asq += (a->e[i])*(a->e[i]);
  asq = (double) 1./sqrt((double)asq);
  for(i=0; i<4; ++i) a->e[i] *= asq;
} 
