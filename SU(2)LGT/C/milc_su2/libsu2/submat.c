/********************  submat.c (in su2.a)  *****************************
*									*
*  Subtract two SU2 matrices 						*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void sub_su2_matrix(su2_matrix *a, su2_matrix *b, su2_matrix *c) {
register int i;
    for(i=0;i<4;i++) c->e[i] = a->e[i] - b->e[i];
}
