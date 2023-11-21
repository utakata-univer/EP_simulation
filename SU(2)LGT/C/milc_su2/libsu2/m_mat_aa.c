/*******************  m_mat_aa.c  (in su2.a) ****************************
*									*
* void mult_su2_aa(su2_matrix *a, su2_matrix *b, su2_matrix *c)		*
* matrix multiply, with both adjoints 						*
* C  <-  A^dagger * B^dagger						*
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_su2_aa(su2_matrix *a, su2_matrix *b, su2_matrix *c) {
  c->e[0] = a->e[0]*b->e[0] - 
            (a->e[1]*b->e[1] + a->e[2]*b->e[2] + a->e[3]*b->e[3]);
  c->e[1] = -a->e[0]*b->e[1] - a->e[1]*b->e[0] + 
            -a->e[2]*b->e[3] + a->e[3]*b->e[2];
  c->e[2] = -a->e[0]*b->e[2] + a->e[1]*b->e[3] - 
            a->e[2]*b->e[0] - a->e[3]*b->e[1];
  c->e[3] = -a->e[0]*b->e[3] - a->e[1]*b->e[2] + 
            a->e[2]*b->e[1] - a->e[3]*b->e[0];
}

