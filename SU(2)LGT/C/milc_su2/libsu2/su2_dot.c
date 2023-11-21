/******************  su2_dot.c  (in su2.a) ******************************
*                                                                       *
* complex su3_dot(su2_vector *a, su2_vector *b);                        *
* return dot product of two su2_vectors                                 *
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

complex su2_dot(su2_vector *a, su2_vector *b) {
   register complex temp1, temp2;
   CMULJ_(a->c[0],b->c[0],temp1);
   CMULJ_(a->c[1],b->c[1],temp2);
   CSUM(temp1,temp2);
   return(temp1);
}
