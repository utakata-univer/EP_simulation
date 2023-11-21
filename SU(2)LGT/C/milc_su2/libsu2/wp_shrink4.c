/*****************  wp_shrink4.c  (in su2.a) ****************************
*									*
* Shrink a wilson vector in four directions, producing four		*
*  half_su2_wilson_vectors.						*
* void wp_shrink_4dir(su2_wilson_vector *a, 
                      half_su2_wilson_vector *b1,
		      half_su2_wilson_vector *b2,
		      half_su2_wilson_vector *b3,
		      half_su2_wilson_vector *b4, int sign)		* 
* B1 <- (1 +- gamma_x)A,, projection					*
*  argument "sign" is sign of gamma matrix.				*
*  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	*
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

void wp_shrink_4dir(su2_wilson_vector *a, 
                      su2_half_wilson_vector *b1,
		      su2_half_wilson_vector *b2,
		      su2_half_wilson_vector *b3,
		      su2_half_wilson_vector *b4, int sign) {
    wp_shrink( a,b1,XUP,sign);
    wp_shrink( a,b2,YUP,sign);
    wp_shrink( a,b3,ZUP,sign);
    wp_shrink( a,b4,TUP,sign);
}
