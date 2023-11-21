/*****************  su3_proj.c  (in su3.a) ******************************
*									*
* void su3_projector(a,b,c) su3_vector *a,*b; su3_matrix *c;		*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "complex.h"
#include "su3.h"

#ifndef FAST
void su3_projector(a,b,c) su3_vector *a,*b; su3_matrix *c; {
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CMUL_J( a->c[i], b->c[j], c->e[i][j] );
    }
}

#else
void su3_projector(a,b,c) su3_vector *a,*b; su3_matrix *c; {
register int i,j;
register float tmp,tmp2;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	tmp2 = a->c[i].real * b->c[j].real;
	tmp = a->c[i].imag * b->c[j].imag;
	c->e[i][j].real = tmp + tmp2;
	tmp2 = a->c[i].real * b->c[j].imag;
	tmp = a->c[i].imag * b->c[j].real;
	c->e[i][j].imag = tmp - tmp2;
    }
}
#endif /* end ifdef FAST */
