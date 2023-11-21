/*******************  cs_m_a_mat.c (in su2.a)  **************************
*									*
*  c_scalar_mult_add_su2mat(su2_matrix *a, su2_matrix *b,               * 
*				  complex phase, su2_matrix *c):        *
*  multiply an su2 matrix by a complex scalar and add it to another     *
*  matrix:   c <- a + phase*b                                           *
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void c_scalar_mult_add_su2_matrix(su2_matrix *a, su2_matrix *b, 
				              complex phase, su2_matrix *c) {
  c->e[0] = a->e[0] + phase.real*b->e[0] - phase.imag*b->e[3];
  c->e[1] = a->e[1] + phase.imag*b->e[0] + phase.real*b->e[3];
  c->e[2] = a->e[2] + phase.real*b->e[2] - phase.imag*b->e[1];
  c->e[3] = a->e[3] + phase.real*b->e[1] + phase.imag*b->e[2];
}
