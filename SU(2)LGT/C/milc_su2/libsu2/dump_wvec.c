/****************  dump_wvec.c  (in su2.a) ***********************
*									*
*  void dump_su2_wvec(su2_wilson_vector *v )				*
*  Print out a Wilson vector 						*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void dump_su2_wvec(su2_wilson_vector *v) {
   register int i,j;
   for(i=0;i<4;i++){
      for(j=0;j<2;j++) printf("(%.2e,%.2e)\t",
	                       v->d[i].c[j].real,v->d[i].c[j].imag);
      printf("\n");
   }
   printf("\n");
}
