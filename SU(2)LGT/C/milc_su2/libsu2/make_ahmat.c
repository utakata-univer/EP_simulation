/*****************  make_ahmat.c  (in su2.a) ****************************
*									*
* void make_anti_hermitian( su2_matrix *m2, anti_hermitmat *ah2)	*
* take the traceless and anti_hermitian part of an su2 matrix 		*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void make_anti_hermitian(su2_matrix *m2, su2_anti_hermitmat *ah2) {

  *ah2 = *m2;
  ah2->e[0] = 0;  
}

