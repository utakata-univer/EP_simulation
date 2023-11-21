/********************  copy_wvec.c  (in su2.a) ********************
*
*void copy_su2_wvec(su2_wilson_vector *src, su2_wilson_vector *dest)
*  copy a Wilson vector
* dest  <-  src
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void copy_su2_wvec(su2_wilson_vector *src, su2_wilson_vector *dest) {
    *dest = *src;	/* hardly worth a function */
}
