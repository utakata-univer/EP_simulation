/*****************  su2_anti_proj.c  (in su2.a) *************************
*									*
* void su2_antiherm_projector(su2_vector *u, su2_vector *v, su2_matrix *m);
* M  <- outer product of U and V, projected onto traceless antihermitian M
*  M_ij = U_i * V_adjoint_j |_projected					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void su2_antiherm_projector(su2_vector *u, su2_vector *v, su2_matrix *m) {

   m->e[0] = 0;
   m->e[1] = 0.5*(u->c[0].imag*v->c[1].real - u->c[0].real*v->c[1].imag
                + u->c[1].imag*v->c[0].real - u->c[1].real*v->c[0].imag);
   m->e[2] = 0.5*(u->c[0].real*v->c[1].real + u->c[0].imag*v->c[1].imag
                - u->c[1].real*v->c[0].real - u->c[1].imag*v->c[0].imag);
   m->e[3] = 0.5*(u->c[0].imag*v->c[0].real - u->c[0].real*v->c[0].imag
                - u->c[1].imag*v->c[1].real + u->c[1].real*v->c[1].imag);
}


