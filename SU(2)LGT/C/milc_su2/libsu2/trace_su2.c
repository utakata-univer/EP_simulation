/*******************  trace_su2.c  (in su3.a) ***************************
*									*
* float trace_su2(a) su2_matrix *a;					*
* return (real) trace of an SU2 matrix 				*
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

/* trace of an SU2 matrix */
float trace_su2(su2_matrix *a) {
    return(2.*a->e[0]);
}
