/************* mb_gamma_r.c  (in su2.a) **************************/
/* 
  Multiply a Wilson matrix by a gamma matrix acting on the column index
  (This is the second index, or equivalently, multiplication on the right)

  usage:  mult_su2_by_gamma_right(su2_wilson_matrix *src, 
                              su2_wilson_matrix *dest, int dir)

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

void mult_su2_by_gamma_right(su2_wilson_matrix *src, 
                              su2_wilson_matrix *dest, int dir) {
register int i; /*color*/
register int c1,s1;	/* row indices, color and spin */

  switch(dir){
    case XUP:
	for(i=0;i<2;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<2;c1++){
	    TIMESMINUSI( src->d[s1].c[c1].d[3].c[i],
		dest->d[s1].c[c1].d[0].c[i] );
	    TIMESMINUSI( src->d[s1].c[c1].d[2].c[i],
		dest->d[s1].c[c1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].c[c1].d[1].c[i],
		dest->d[s1].c[c1].d[2].c[i] );
	    TIMESPLUSI(  src->d[s1].c[c1].d[0].c[i],
		dest->d[s1].c[c1].d[3].c[i] );
	}
	break;
    case YUP:
	for(i=0;i<2;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<2;c1++){
	    TIMESMINUSONE( src->d[s1].c[c1].d[3].c[i],
		dest->d[s1].c[c1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].c[c1].d[2].c[i],
		dest->d[s1].c[c1].d[1].c[i] );
	    TIMESPLUSONE(  src->d[s1].c[c1].d[1].c[i],
		dest->d[s1].c[c1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].c[c1].d[0].c[i],
		dest->d[s1].c[c1].d[3].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<2;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<2;c1++){
	    TIMESMINUSI( src->d[s1].c[c1].d[2].c[i],
		dest->d[s1].c[c1].d[0].c[i] );
	    TIMESPLUSI(  src->d[s1].c[c1].d[3].c[i],
		dest->d[s1].c[c1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].c[c1].d[0].c[i],
		dest->d[s1].c[c1].d[2].c[i] );
	    TIMESMINUSI( src->d[s1].c[c1].d[1].c[i],
		dest->d[s1].c[c1].d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<2;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<2;c1++){
	    TIMESPLUSONE( src->d[s1].c[c1].d[2].c[i],
		dest->d[s1].c[c1].d[0].c[i] );
	    TIMESPLUSONE( src->d[s1].c[c1].d[3].c[i],
		dest->d[s1].c[c1].d[1].c[i] );
	    TIMESPLUSONE( src->d[s1].c[c1].d[0].c[i],
		dest->d[s1].c[c1].d[2].c[i] );
	    TIMESPLUSONE( src->d[s1].c[c1].d[1].c[i],
		dest->d[s1].c[c1].d[3].c[i] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<2;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<2;c1++){
	    TIMESPLUSONE(  src->d[s1].c[c1].d[0].c[i],
		dest->d[s1].c[c1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].c[c1].d[1].c[i],
		dest->d[s1].c[c1].d[1].c[i] );
	    TIMESMINUSONE( src->d[s1].c[c1].d[2].c[i],
		dest->d[s1].c[c1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].c[c1].d[3].c[i],
		dest->d[s1].c[c1].d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
  }
}

