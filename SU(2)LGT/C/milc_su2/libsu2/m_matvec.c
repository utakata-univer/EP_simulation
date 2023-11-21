/****************  m_matvec.c  (in su2.a) *******************************
*									*
* void mult_su2_mat_vec(su2_matrix *a, su2_vector *b, su2_vector *c)	*
* matrix times vector multiply, no adjoints 				*
*  C  <-  A*B								*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_su2_mat_vec(su2_matrix *a, su2_vector *b, su2_vector *c) {

   c->c[0].real =  a->e[0]*b->c[0].real
                  - a->e[3]*b->c[0].imag
                  + a->e[2]*b->c[1].real
                  - a->e[1]*b->c[1].imag;

   c->c[0].imag = a->e[3]*b->c[0].real
                  + a->e[0]*b->c[0].imag
                  + a->e[2]*b->c[1].imag
                  + a->e[1]*b->c[1].real;

   c->c[1].real = - a->e[2]*b->c[0].real
                  - a->e[1]*b->c[0].imag
                  + a->e[0]*b->c[1].real
                  + a->e[3]*b->c[1].imag;

   c->c[1].imag = - a->e[2]*b->c[0].imag
                  + a->e[1]*b->c[0].real
                  + a->e[0]*b->c[1].imag
                  - a->e[3]*b->c[1].real;
}



