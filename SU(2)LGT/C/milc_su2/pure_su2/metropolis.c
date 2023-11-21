/************************** metropolis.c *******************************/
/* Metropolis updating for SU2 pure gauge */
/* MIMD version 4 */
/* J. Hetrick and D. Toussaint April 1995 */


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

/* Generic definitions - could be useful elsewhere */
#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/


void monte(int NumStp) {
  int NumTrj,Nhit;
  int parity;
  float scale;		/* limits size of change matrix */
  su2_matrix change;	/* proposed change in link */
  su2_matrix newlink;	/* change * oldlink */
  int dir, i;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  float oldaction,newaction;

  accept = reject = 0;
  scale = 1.0;
  for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++) {

    for(parity=ODD;parity<=EVEN;parity++) {
      FORALLUPDIR(dir) {
	/* compute the gauge force */
	dsdu_qhb(dir,parity); 
	/* now for the Metropolis updating */
	FORSOMEPARITY(i,st,parity) {
	  /* generate random SU(2) matrix */
	  /* scale < 2/sqrt(3), so vector magnitude < 1 */
	  change.e[1] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[2] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[3] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[0] = sqrt(1.0 - change.e[1]*change.e[1] -
		     change.e[2]*change.e[2] - change.e[3]*change.e[3]);
	  mult_su2_nn( &change, &(st->link[dir]), &newlink );
	  
	  /* compute old action and new action */
	  oldaction=(0.5*beta)*realtrace_su2( &(st->link[dir]), &(st->staple) );
	  newaction=(0.5*beta)*realtrace_su2( &newlink, &(st->staple) );

	  /* accept or reject */
	  if( newaction > oldaction ){
	    st->link[dir]=newlink;
	    accept++;
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	      st->link[dir]=newlink;
	      accept++;
	    }
	    else{
	      reject++;
	   }
	  }
      
	} /*   st */

      } /*  direction */
    }
  } /* parity and NumTrj */
  /* diagnostics: */
/*  printf("monte: accept = %d, reject = %d\n",accept,reject); */
  
} /* monte */




