/*****************  su3_projector_w.c  (in su3.a) ******************************
*									*
* void su3_projector_w(a,b,c) wilson_vector *a,*b; su3_matrix *c;	*
* C  <- sum over spins of outer product of A.d[i] and B.d[i]		*
*  C_ij = sum( A_i * B_adjoint_j )					*
*/
#include "complex.h"
#include "su3.h"

#ifndef FAST
void su3_projector_w(a,b,c) wilson_vector *a,*b; su3_matrix *c; {
register int i,j,k;
register complex cc;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	c->e[i][j] = cmplx(0.0,0.0);
	for(k=0;k<4;k++){
	    CMUL_J( a->d[k].c[i], b->d[k].c[j], cc ); CSUM( c->e[i][j], cc );
	}
    }
}

#else
void su3_projector_w(a,b,c) wilson_vector *a,*b; su3_matrix *c; {
register int i,j,k;
register float tmp_r,tmp_i,tmp2;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	tmp_r = tmp_i = 0.0;
	for(k=0;k<4;k++){
	    tmp2 = a->d[k].c[i].real * b->d[k].c[j].real; tmp_r = tmp_r + tmp2;
	    tmp2 = a->d[k].c[i].imag * b->d[k].c[j].imag; tmp_r = tmp_r + tmp2;
	    tmp2 = a->d[k].c[i].imag * b->d[k].c[j].real; tmp_i = tmp_i + tmp2;
	    tmp2 = a->d[k].c[i].real * b->d[k].c[j].imag; tmp_i = tmp_i - tmp2;
	}

	c->e[i][j].real = tmp_r;
	c->e[i][j].imag = tmp_i;
    }
}
#endif /* end ifdef FAST */
