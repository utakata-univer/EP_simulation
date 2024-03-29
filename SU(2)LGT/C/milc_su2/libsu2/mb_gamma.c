/************* mb_gamma.c  (in su2.a) **************************/
/* 
  Multiply a Wilson vector by a gamma matrix

void mult_su2_by_gamma(su2_wilson_vector *src, su2_wilson_vector *dest,
		       int dir );

  usage:  mult_su2_by_gamma(...)

          dir = XUP, YUP, ZUP, TUP or GAMMAFIVE

 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 
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

void mult_su2_by_gamma(su2_wilson_vector *src, su2_wilson_vector *dest,
		       int dir ) {
  register int i; /*color*/

  switch(dir){
    case XUP:
	for(i=0;i<2;i++){
	    TIMESPLUSI(  src->d[3].c[i], dest->d[0].c[i] );
	    TIMESPLUSI(  src->d[2].c[i], dest->d[1].c[i] );
	    TIMESMINUSI( src->d[1].c[i], dest->d[2].c[i] );
	    TIMESMINUSI( src->d[0].c[i], dest->d[3].c[i] );
	}
	break;
    case YUP:
	for(i=0;i<2;i++){
	    TIMESMINUSONE( src->d[3].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE(  src->d[2].c[i], dest->d[1].c[i] );
	    TIMESPLUSONE(  src->d[1].c[i], dest->d[2].c[i] );
	    TIMESMINUSONE( src->d[0].c[i], dest->d[3].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<2;i++){
	    TIMESPLUSI(  src->d[2].c[i], dest->d[0].c[i] );
	    TIMESMINUSI( src->d[3].c[i], dest->d[1].c[i] );
	    TIMESMINUSI( src->d[0].c[i], dest->d[2].c[i] );
	    TIMESPLUSI(  src->d[1].c[i], dest->d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<2;i++){
	    TIMESPLUSONE( src->d[2].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE( src->d[3].c[i], dest->d[1].c[i] );
	    TIMESPLUSONE( src->d[0].c[i], dest->d[2].c[i] );
	    TIMESPLUSONE( src->d[1].c[i], dest->d[3].c[i] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<2;i++){
	    TIMESPLUSONE(  src->d[0].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE(  src->d[1].c[i], dest->d[1].c[i] );
	    TIMESMINUSONE( src->d[2].c[i], dest->d[2].c[i] );
	    TIMESMINUSONE( src->d[3].c[i], dest->d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA()\n");
  }
}

