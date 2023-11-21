/*********************** reunitarize.c ***************************/
/* MIMD version 3 */

/* reunitarize the link matrices */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>
#define TOLERANCE (0.0001)
/*#define UNIDEBUG */
float check_su2();

void check_unitarity() {
  register int i,dir;
  int ii,jj;
  register site *s;
  register su2_matrix *mat;
  float deviation,max_deviation;
  double av_deviation;
  union {
    float fval;
    int ival;
  } xxx;

  max_deviation=av_deviation=0;
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; dir++ ){
      mat = (su2_matrix *)&(s->link[dir]);
      deviation=check_su2( mat );
      if (deviation>TOLERANCE){
	printf("Unitarity problem on node %d, site %d, dir %d, deviation=%f\n",
	       mynode(),i,dir,deviation);
	printf("SU2 matrix:\n");
	dump_su2_mat(mat);

	if(max_deviation<deviation) max_deviation=deviation;
	av_deviation += deviation*deviation;
      }
    }
  }
  av_deviation = sqrt(av_deviation/(4*i));

#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %g, avrg %g\n",
	 mynode(), max_deviation, av_deviation);
#endif

  if(max_deviation> TOLERANCE) 
    printf("Unitarity problem on node %d, maximum deviation=%f\n",
	   mynode(),max_deviation);
}  /*check_unitarity() */


float check_su2(su2_matrix *c) { return(fabs((double)det_su2(c) - 1.)); }
