/********************  msq_wvec.c  (in su2.a) ********************
*
* float msq_su2_wvec(su2_wilson_vector *vec)
*  squared magnitude of a Wilson vector
* 
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float magsq_su2_wvec(su2_wilson_vector *vec) {
  register int i;
  register float sum;
  for(i=0,sum=0.0;i<4;i++) sum += magsq_su2vec( &(vec->d[i]) );
  return(sum);
}
