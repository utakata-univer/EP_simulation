/*****************  cmp_ahmat.c  (in su2.a) ****************************
*									*
* void compress_anti_hermitian( su2_matrix *m2, su2_anti_hermitmat *ah2)	*
* does nothing in su2; copies m2 -> ah2
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void compress_anti_hermitian(su2_matrix *m2, su2_anti_hermitmat *ah2) {

  *ah2 = *m2;
}

