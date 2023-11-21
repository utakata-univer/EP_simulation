/*****************  su2mat_copy.c  (in su2.a) ***************************
*									*
* void su2mat_copy(su2_matrix *a, su2_matrix *b);			*
* Copy an su2 matrix:  B <- A   					*
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void su2mat_copy(su2_matrix *a, su2_matrix *b) {
  register int i;
 
  *b = *a;
}
