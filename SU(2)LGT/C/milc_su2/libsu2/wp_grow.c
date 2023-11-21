/***************** wp_grow.c  (in su2.a) **************************/
/* 
  Expand the "Wilson projection" of a Wilson fermion vector.
  (1 +- gamma_j) is a projection operator, and we are given a
  half_wilson_vector which contains the two components of a Wilson
  vector projected out.  This routine reexpands it to a four component
  object.

  void su2_wp_grow(su2_half_wilson_vector *src, su2_wilson_vector *dest,
		 int dir, int sign);

  usage:  su2_wp_grow(...)

	If dir is one of XUP,YUP,ZUP or TUP, the projection is
	along the eigenvectors with eigenvalue +1, which survive
	multiplication by (1+gamma[dir]).
	If dir is one of XDOWN,YDOWN,ZDOWN or TDOWN, the projection is
	along the eigenvectors with eigenvalue -1, which survive
	multiplication by (1-gamma[OPP_DIR(dir)]).
	If sign=MINUS reverse the roles of +1 and -1 - in other words
	use -gamma_dir instead of gamma_dir

  Here my eigenvectors are normalized to 2, so for XYZT directions
  I won't explicitely multiply by 2.  In other words, the matrix of
  eigenvectors is sqrt(2) times a unitary matrix, and in reexpanding
  the vector I will multiply by the adjoint of this matrix.

  For UP directions, hvec.h[0] and hvec.h[2] contain the projections
  along the first and second eigenvectors respectively.
  For DOWN directions, hvec.h[0] and hvec.h[2] contain the projections
  along the third and fourth eigenvectors respectively. This results
  in down directions differing from up directions only in the sign of
  the addition.

  Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to multiplication
   by 1+-gamma_dir

 gamma(XUP) 			eigenvectors	eigenvalue
 	    0  0  0  i		( 1, 0, 0,-i)	+1
            0  0  i  0		( 0, 1,-i, 0)	+1
            0 -i  0  0		( 0, 1, 0,+i)	-1
           -i  0  0  0		( 1, 0,+i, 0)	-1

 gamma(YUP)			eigenvectors	eigenvalue
 	    0  0  0 -1		( 1, 0, 0,-1)	+1
            0  0  1  0		( 0, 1, 1, 0)	+1
            0  1  0  0		( 1, 0, 0, 1)	-1
           -1  0  0  0		( 0, 1,-1, 0)	-1

 gamma(ZUP)			eigenvectors	eigenvalue
 	    0  0  i  0		( 1, 0,-i, 0)	+1
            0  0  0 -i		( 0, 1, 0,+i)	+1
           -i  0  0  0		( 1, 0,+i, 0)	-1
            0  i  0  0		( 0, 1, 0,-i)	-1

 gamma(TUP)			eigenvectors	eigenvalue
 	    0  0  1  0		( 1, 0, 1, 0)	+1
            0  0  0  1		( 0, 1, 0, 1)	+1
            1  0  0  0		( 1, 0,-1, 0)	-1
            0  1  0  0		( 0, 1, 0,-1)	-1

 gamma(FIVE) 			eigenvectors	eigenvalue
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
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

void su2_wp_grow(su2_half_wilson_vector *src, su2_wilson_vector *dest,
		 int dir, int sign) {
  register int i; /*color*/

  if(sign==MINUS)dir=OPP_DIR(dir);	/* two ways to get -gamma_dir ! */
  switch(dir){
    case XUP:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESMINUSI( src->h[0].c[i], dest->d[3].c[i]);
	    TIMESMINUSI( src->h[1].c[i], dest->d[2].c[i]);
	}
	break;
    case XDOWN:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESPLUSI( src->h[0].c[i], dest->d[3].c[i]);
	    TIMESPLUSI( src->h[1].c[i], dest->d[2].c[i]);
	}
	break;
    case YUP:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESMINUSONE( src->h[0].c[i], dest->d[3].c[i] );
	    TIMESPLUSONE(  src->h[1].c[i], dest->d[2].c[i] );
	}
	break;
    case YDOWN:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESPLUSONE(  src->h[0].c[i], dest->d[3].c[i] );
	    TIMESMINUSONE( src->h[1].c[i], dest->d[2].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESMINUSI( src->h[0].c[i], dest->d[2].c[i] );
	    TIMESPLUSI(  src->h[1].c[i], dest->d[3].c[i] );
	}
	break;
    case ZDOWN:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESPLUSI(  src->h[0].c[i], dest->d[2].c[i] );
	    TIMESMINUSI( src->h[1].c[i], dest->d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    dest->d[2].c[i]      = src->h[0].c[i];
	    dest->d[3].c[i]      = src->h[1].c[i];
	}
	break;
    case TDOWN:
	for(i=0;i<2;i++){
	    dest->d[0].c[i]      = src->h[0].c[i];
	    dest->d[1].c[i]      = src->h[1].c[i];
	    TIMESMINUSONE( src->h[0].c[i], dest->d[2].c[i] );
	    TIMESMINUSONE( src->h[1].c[i], dest->d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO WP_GROW()\n");
  }
}

