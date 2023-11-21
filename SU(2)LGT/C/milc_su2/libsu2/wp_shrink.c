/************* wp_shrink.c  (in su2.a) **************************/
/* 

  Compute the "Wilson projection" of a Wilson fermion vector.
  (1 +- gamma_j) is a projection operator, and we want to isolate
  the components of the vector that it keeps.  In other words, keep
  the components of the vector along the eigenvectors of 1+-gamma_j
  with eigenvalue 2, and throw away those with eigenvalue 0.

  void su2_wp_shrink(su2_wilson_vector *src, half_su2_wilson_vector *dest,
		   int dir, int sign);

  usage:  su2_wp_shrink( src, dest, dir, sign)
	wilson_vector *src;
	half_wilson_vector *dest;
	int dir,sign;

	If dir is one of XUP,YUP,ZUP or TUP, take the projections
	along the eigenvectors with eigenvalue +1, which survive
	multiplication by (1+gamma[dir]).
	If dir is one of XDOWN,YDOWN,ZDOWN or TDOWN, take the projections
	along the eigenvectors with eigenvalue -1, which survive
	multiplication by (1-gamma[OPP_DIR(dir)]).
	If sign=MINUS, switch the roles of +1 and -1 (ie use -gamma_dir
	instead of gamma_dir )

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

  Note: su2_wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to multiplication
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

void su2_wp_shrink(su2_wilson_vector *src, su2_half_wilson_vector *dest,
		   int dir, int sign) {
  register int i; /*color*/

  if(sign==MINUS)dir=OPP_DIR(dir);	/* two ways to get -gamma_dir ! */
  switch(dir){
    case XUP:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[3].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[3].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[2].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[2].c[i].real;
	}
	break;
    case XDOWN:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[3].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[3].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[2].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[2].c[i].real;
	}
	break;
    case YUP:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[3].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[3].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[2].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[2].c[i].imag;
	}
	break;
    case YDOWN:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[3].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[3].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[2].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[2].c[i].imag;
	}
	break;
    case ZUP:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[2].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[2].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[3].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[3].c[i].real;
	}
	break;
    case ZDOWN:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[2].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[2].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[3].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[3].c[i].real;
	}
	break;
    case TUP:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[2].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[2].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[3].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[3].c[i].imag;
	}
	break;
    case TDOWN:
	for(i=0;i<2;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[2].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[2].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[3].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[3].c[i].imag;
	}
	break;
    default:
	printf("BAD CALL TO SU2_WP_SHRINK()\n");
  }
}

