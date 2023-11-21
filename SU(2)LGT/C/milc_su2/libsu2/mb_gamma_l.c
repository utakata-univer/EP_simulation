/************* mb_gamma_l.c  (in su2.a) **************************/
/* 
  Multiply a Wilson matrix by a gamma matrix acting on the row index
  (This is the first index, or equivalently, multiplication on the left)

void mult_su2_by_gamma_left(su2_wilson_matrix *src, su2_wilson_matrix *dest,
			    int dir);

  usage:  mult_su2_by_gamma_left(...)

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

void mult_su2_by_gamma_left(su2_wilson_matrix *src, su2_wilson_matrix *dest,
			    int dir) {
  register int i; /*color*/
  register int c2,s2;	/* column indices, color and spin */

  switch(dir){
    case XUP:
	for(i=0;i<2;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<2;c2++){
	    TIMESPLUSI(  src->d[3].c[i].d[s2].c[c2],
		dest->d[0].c[i].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[2].c[i].d[s2].c[c2],
		dest->d[1].c[i].d[s2].c[c2] );
	    TIMESMINUSI( src->d[1].c[i].d[s2].c[c2],
		dest->d[2].c[i].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].c[i].d[s2].c[c2],
		dest->d[3].c[i].d[s2].c[c2] );
	}
	break;
    case YUP:
	for(i=0;i<2;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<2;c2++){
	    TIMESMINUSONE( src->d[3].c[i].d[s2].c[c2],
		dest->d[0].c[i].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[2].c[i].d[s2].c[c2],
		dest->d[1].c[i].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].c[i].d[s2].c[c2],
		dest->d[2].c[i].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[0].c[i].d[s2].c[c2],
		dest->d[3].c[i].d[s2].c[c2] );
	}
	break;
    case ZUP:
	for(i=0;i<2;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<2;c2++){
	    TIMESPLUSI(  src->d[2].c[i].d[s2].c[c2],
		dest->d[0].c[i].d[s2].c[c2] );
	    TIMESMINUSI( src->d[3].c[i].d[s2].c[c2],
		dest->d[1].c[i].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].c[i].d[s2].c[c2],
		dest->d[2].c[i].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[1].c[i].d[s2].c[c2],
		dest->d[3].c[i].d[s2].c[c2] );
	}
	break;
    case TUP:
	for(i=0;i<2;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<2;c2++){
	    TIMESPLUSONE( src->d[2].c[i].d[s2].c[c2],
		dest->d[0].c[i].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[3].c[i].d[s2].c[c2],
		dest->d[1].c[i].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[0].c[i].d[s2].c[c2],
		dest->d[2].c[i].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[1].c[i].d[s2].c[c2],
		dest->d[3].c[i].d[s2].c[c2] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<2;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<2;c2++){
	    TIMESPLUSONE(  src->d[0].c[i].d[s2].c[c2],
		dest->d[0].c[i].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].c[i].d[s2].c[c2],
		dest->d[1].c[i].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[2].c[i].d[s2].c[c2],
		dest->d[2].c[i].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[3].c[i].d[s2].c[c2],
		dest->d[3].c[i].d[s2].c[c2] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_LEFT()\n");
  }
}

