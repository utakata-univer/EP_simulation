/******************  dump_su2_mat.c  (in su2.a) *************************
*									*
*  void dump_su2_mat( su2_matrix *mat )					*
*  print out a 2x2 SU(2) matrix: a0 a1 a2 a3					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void dump_su2_mat(su2_matrix *m) {
   printf("%f\t%f\t%f\t%f\n",m->e[0],m->e[1],m->e[2],m->e[3]);
}
