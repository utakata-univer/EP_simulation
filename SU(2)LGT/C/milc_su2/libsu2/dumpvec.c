/*******************  dumpvec.c  (in su2.a) *****************************
*									*
*  void dump_su2_vec( su2_vector *vec )					*
*  print out a 2 element complex vector					*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void dump_su2_vec(su2_vector *v) {
   int j;
   for(j=0;j<2;j++) printf("(%.2e,%.2e)\t",
			   v->c[j].real,v->c[j].imag);
   printf("\n");
}
