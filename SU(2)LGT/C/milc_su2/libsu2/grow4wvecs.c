/*****************  grow4wvecs.c  (in su2.a) ****************************
*									*
*  If sum=0,								*
*  Grow and add four su2_wilson_vectors					*
*  If sum=1,								*
*  Grow and sum four su2_wilson_vectors to another su2_wilson_vector	*
* 
* void grow_add_four_su2_wvecs(su2_wilson_vector *a,
* 			     su2_half_wilson_vector *b1,
* 			     su2_half_wilson_vector *b2,
* 			     su2_half_wilson_vector *b3,
* 			     su2_half_wilson_vector *b4, int sign, int sum);
* 
/* * A  <-  B1 + B2 + B3 + B4   or					*
* A  <-  A + B1 + B2 + B3 + B4						*
* B1 is expanded using gamma_x, B2 using gamma_y, etc. 			*
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

/* Directions, and a macro to give the opposite direction */
/*  These must go from 0 to 7 because they will be used to index an
    array. */
/* Also define NDIRS = number of directions */
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define OPP_DIR(dir)	(7-(dir))	/* Opposite direction */
#define NDIRS 8				/* number of directions */

void grow_add_four_su2_wvecs(su2_wilson_vector *a,
			     su2_half_wilson_vector *b1,
			     su2_half_wilson_vector *b2,
			     su2_half_wilson_vector *b3,
			     su2_half_wilson_vector *b4, int sign, int sum) {
    if(sum==0)su2_wp_grow( b1,a,XUP,sign);
    else su2_wp_grow_add( b1,a,XUP,sign);
    su2_wp_grow_add( b2,a,YUP,sign);
    su2_wp_grow_add( b3,a,ZUP,sign);
    su2_wp_grow_add( b4,a,TUP,sign);
}
