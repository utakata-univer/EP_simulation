/***************  m_mat_choose.c  (in su2.a) ****************************
*									*
* void mult_su2_choose(su2_matrix *a, int flaga,                        *
*                      su2_matrix *b, int flagb, su2_matrix *c)   	*
* matrix multiply, adjoints set by flags                                *
* C  <-  A^(flaga) * B^(flagb)						*
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void mult_su2_choose(su2_matrix *a, int flaga, su2_matrix *b, 
		 int flagb, su2_matrix *c) {
  c->e[0] = a->e[0]*b->e[0] - flaga*flagb*(a->e[1]*b->e[1] + 
					 a->e[2]*b->e[2] + a->e[3]*b->e[3]);
  c->e[1] = flagb*a->e[0]*b->e[1] + flaga*a->e[1]*b->e[0] + 
            flaga*flagb*(-a->e[2]*b->e[3] + a->e[3]*b->e[2]);
  c->e[2] = flagb*a->e[0]*b->e[2] + flaga*flagb*a->e[1]*b->e[3] + 
            flaga*a->e[2]*b->e[0] - flaga*flagb*a->e[3]*b->e[1];
  c->e[3] = flagb*a->e[0]*b->e[3] - flaga*flagb*a->e[1]*b->e[2] + 
            flaga*flagb*a->e[2]*b->e[1] + flaga*a->e[3]*b->e[0];
}

