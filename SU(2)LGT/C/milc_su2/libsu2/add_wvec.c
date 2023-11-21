/********************  add_wvec.c  (in su2.a) ********************
*
* void add_su2_wilson_vector(su2_wilson_vector *src1, 
*                 su2_wilson_vector *src2, su2_wilson_vector *dest)
*  add two Wilson vectors
*  dest  <-  src1 + src2
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void add_su2_wilson_vector(su2_wilson_vector *src1, su2_wilson_vector *src2,
                                               su2_wilson_vector *dest) {
  register int i;
  for(i=0;i<4;i++) add_su2_vector( &(src1->d[i]), 
                                   &(src2->d[i]), &(dest->d[i]));
}
