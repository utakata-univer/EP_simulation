/******************  det_su2.c  (in su2.a) ******************************
*									*
* float det_su2( su2_matrix *a )					*
* determinant of an su2_matrix: if a is normalized, det = 1	        *
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float det_su2(su2_matrix *a) {
   return(a->e[0]*a->e[0] + a->e[1]*a->e[1] + a->e[2]*a->e[2] + 
	  a->e[3]*a->e[3]);
}
