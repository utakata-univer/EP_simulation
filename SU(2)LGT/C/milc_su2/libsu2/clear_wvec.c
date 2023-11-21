/********************  clear_wvec.c  (in su2.a) ********************
*
*void clear_su2_wilson_vector(su2_wilson_vector *dest)
*  clear a Wilson vector
*  dest  <-  zero_vector
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void clear_wvec(su2_wilson_vector *dest) {
   register int i,j;
   for(i=0;i<4;i++)for(j=0;j<2;j++){
      dest->d[i].c[j].real = dest->d[i].c[j].imag = 0.0;
   }
}
