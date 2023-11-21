/*****************  sub4vecs.c  (in su2.a) ******************************
*									*
*  Subtract four su2_vectors from an su2_vector				*
* void sub_four_su2_vecs(su2_vector *a, su2_vector *b1, su2_vector *b2, *
*                                  su2_vector *b3, su2_vector *b4)	*
* A  <-  A - B1 - B2 - B3 - B4						*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void sub_four_su2_vecs(su2_vector *a, su2_vector *b1, su2_vector *b2, 
		       su2_vector *b3, su2_vector *b4) {
  register int i;
  for(i=0;i<2;i++) {
     CSUB( a->c[i], b1->c[i], a->c[i] );
     CSUB( a->c[i], b2->c[i], a->c[i] );
     CSUB( a->c[i], b3->c[i], a->c[i] );
     CSUB( a->c[i], b4->c[i], a->c[i] );
  }
}
