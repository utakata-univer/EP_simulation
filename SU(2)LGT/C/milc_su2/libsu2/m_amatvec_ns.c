/****************  m_amatvec_ns.c  (in su2.a) *******************************
*				  					   *
* void mult_adj_su2_mat_vec_nsum(su2_matrix *a, su2_vector *b, su2_vector *c)
* adj su2_matrix times su2_vector multiply and subtract from another       *
*  C  <- C - A_adjoint*B						   *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_adj_su2_mat_vec_nsum(su2_matrix *a, su2_vector *b, 
                                              su2_vector *c) {

   c->c[0].real -=  a->e[0]*b->c[0].real
                  + a->e[3]*b->c[0].imag
                  - a->e[2]*b->c[1].real
                  + a->e[1]*b->c[1].imag;

   c->c[0].imag -= -a->e[3]*b->c[0].real
                  + a->e[0]*b->c[0].imag
                  - a->e[2]*b->c[1].imag
                  - a->e[1]*b->c[1].real;

   c->c[1].real -=  a->e[2]*b->c[0].real
                  + a->e[1]*b->c[0].imag
                  + a->e[0]*b->c[1].real
                  - a->e[3]*b->c[1].imag;

   c->c[1].imag -=  a->e[2]*b->c[0].imag
                  - a->e[1]*b->c[0].real
                  + a->e[0]*b->c[1].imag
                  + a->e[3]*b->c[1].real;
}



