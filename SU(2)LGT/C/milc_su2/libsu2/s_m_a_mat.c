/********************  s_m_a_mat.c (in su2.a)  **************************
*									*
* void scalar_mult_add_su2_matrix(a,b,s,c) su2_matrix *a,*b,*c; float s *
* C <- A + s*B                                                          *
*/

#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void scalar_mult_add_su2_matrix(su2_matrix *a, su2_matrix *b, float s, 
		    su2_matrix *c) {
register int i;
    for(i=0;i<4;i++) c->e[i] = a->e[i] + s*b->e[i];
}
