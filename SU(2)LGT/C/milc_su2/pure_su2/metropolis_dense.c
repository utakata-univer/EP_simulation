
/************************** metropolis_dense.c *******************************/
/* Metropolis updating for SU2 pure gauge */
/* MIMD version 4 */
/* J. Hetrick and D. Toussaint April 1995 */
/* update with "almost quenched" approximation for nonzero density 6/1/95 */
/* monte_space() does spatial links, monte_time() does temporal */


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

/* Generic definitions - could be useful elsewhere */
#define FORSPACEUPDIR(dir) for(dir=XUP; dir<=ZUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/

void ploop_less_slice(int time,int parity);

test_ploop(){
int i,time; site *s; su2_matrix tmat;
  FORALLSITES(i,s){
    s->link[TUP] =  make_su2_matrix( (float)(s->x),0.0,0.0,0.0 );
  }
  for(time=0;time<nt;time++){
    FORALLSITES(i,s){
      s->ploop = make_su2_matrix(99.9,0.0,0.0,0.0);
    }
    ploop_less_slice(time,EVEN);
    FORALLSITES(i,s)if(s->t==time){
      printf("site %d, coords %d %d %d %d\n",i,s->x,s->y,s->z,s->t);
      dump_su2_mat( &(s->ploop) );
    }
  }
}


void monte_space(int NumStp) {
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
      FORSPACEUPDIR(dir) {
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
  /* diagnostics */
  printf("monte_space: accept = %d, reject = %d\n",accept,reject);
} /* monte_space */





/* time direction includes det(Polyakov+C) in its action */
void monte_time(int NumStp) {
  int NumTrj,Nhit;
  int time,parity;
  float scale;		/* limits size of change matrix */
  float C;		/* constant to add to ploop matrix */
  su2_matrix change;	/* proposed change in link */
  su2_matrix newlink;	/* change * oldlink */
  su2_matrix pc;	/* "almost Polyakov loop" + constant */
  su2_matrix tmat;	/* scratch */
  int i;
  register site *st;
  int accept, reject;	/* number of accepts and rejects */
  float oldaction,newaction;

  accept = reject = 0;
  scale = 1.0;
  C = 1.0; /*FIXXX ME*/
  for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++) {

    for(parity=ODD;parity<=EVEN;parity++) {
      /* compute the gauge force */
      dsdu_qhb(TUP,parity); 
      /* update time slices separately, because they are linked by P. loops */
      for(time=0;time<nt;time++){
	ploop_less_slice(time,parity);
	/* now for the Metropolis updating */
	FORSOMEPARITY(i,st,parity)if(st->t==time) {
	  /* generate random SU(2) matrix */
	  /* scale < 2/sqrt(3), so vector magnitude < 1 */
	  change.e[1] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[2] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[3] = scale *( myrand(&(st->site_prn)) - 0.5 );
	  change.e[0] = sqrt(1.0 - change.e[1]*change.e[1] -
		     change.e[2]*change.e[2] - change.e[3]*change.e[3]);
	  mult_su2_nn( &change, &(st->link[TUP]), &newlink );
	  
	  /* compute old action and new action */
	  mult_su2_nn( &(st->link[TUP]), &(st->ploop), &tmat );
	  tmat.e[0] += C;  /* add C*identity matrix to matrix (SU2 trick)*/
	  oldaction=(0.5*beta)*realtrace_su2( &(st->link[TUP]), &(st->staple) );
	  oldaction += log( det_su2( &tmat ) );
/*
printf("OLD: site = %d = %d %d %d %d    %e\n",
i,st->x,st->y,st->z,st->t,log( det_su2( &tmat ) ) );
dump_su2_mat( &(st->link[TUP]) );
dump_su2_mat( &(st->ploop) );
dump_su2_mat( &tmat );
*/

	  mult_su2_nn( &newlink, &(st->ploop), &tmat );
	  tmat.e[0] += C;  /* add C*identity matrix to matrix (SU2 trick)*/
	  newaction=(0.5*beta)*realtrace_su2( &newlink, &(st->staple) );
	  newaction += log( det_su2( &tmat ) );
/*
printf("NEW: site = %d = %d %d %d %d    %e\n",
i,st->x,st->y,st->z,st->t,log( det_su2( &tmat ) ) );
dump_su2_mat( &newlink );
dump_su2_mat( &(st->ploop) );
dump_su2_mat( &tmat );
printf("\n");
*/
	  /* accept or reject */
	  if( newaction > oldaction ){
	    st->link[TUP]=newlink;
	    accept++;
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	      st->link[TUP]=newlink;
	      accept++;
	    }
	    else{
	      reject++;
	   }
	  }
      
	} /*   st */

      } /* time slice */
    } /* parity */
  } /* NumTrj */
  /* diagnostics */
  printf("monte_time: accept = %d, reject = %d\n",accept,reject);
} /* monte_time */
	

/* Calculate product of timelike links for all time slices except "time",
   for sites where parity at "time" is "parity".
   Put the result in the matrix "ploop"
*/
void ploop_less_slice(int time,int parity){
  int l_time;	/* time at which we are multiplying */
  int c_time;   /* l_time%nt - use this one in accessing sites */
  int l_parity; 
	/* parity of sites at l_time where parity at "time" is "parity" */
  register int i;
  register site *s;
  msg_tag *tag;

  if(parity==EVENANDODD){
    if(this_node==0)printf("Bad parity in ploop_less_slice()\n");
    terminate(0);
  }

  FORSOMEPARITY( i,s,OPP_PAR(parity) )if(s->t==(time-1+nt)%nt){
    s->ploop = s->link[TUP];
  }

  for( l_time = time + nt-2; l_time > time; l_time--){
    c_time = l_time >= nt ? l_time-nt : l_time;
    if( (l_time-time)%2 ==0){l_parity=parity;}
    else {l_parity = OPP_PAR(parity);}

    /* gather current product from slice above */
    tag = start_gather( F_OFFSET(ploop), sizeof(su2_matrix), TUP,
      l_parity, gen_pt[0]);
    wait_gather(tag);
    FORSOMEPARITY( i,s,l_parity )if(s->t==c_time){
      mult_su2_nn( &(s->link[TUP]), (su2_matrix *)gen_pt[0][i],
	&(s->ploop) );
	/* since only on one time slice, don't have to worry if
	   gen_pt points to ploop */
    cleanup_gather(tag);
    } /* end loop over sites at time l_time */
  } /* end loop over l_times */

  /* finish by bringing result to desired time slice */
  tag = start_gather( F_OFFSET(ploop), sizeof(su2_matrix), TUP,
    parity, gen_pt[0]);
  wait_gather(tag);
  FORSOMEPARITY( i,s,parity )if(s->t==time)
    s->ploop = *(su2_matrix *)gen_pt[0][i];
  cleanup_gather(tag);
}
