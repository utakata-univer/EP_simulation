/*****************  su2vec_copy.c  (in su2.a) ***************************
*									*
* void su2vec_copy(su2_vector *a, su2_vector *b);			*
* Copy an su2 vector:  B <- A   					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void su2vec_copy(su2_vector *a, su2_vector *b) { *b = *a; }
