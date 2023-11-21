/****************** ploop2.c ************************************/
/* MIMD version 3 */
/* evaluate the Polyakov loops.  This version uses gathers. */

#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "su3.h"
#include LATDEF
#include <comdefs.h>

complex ploop() {
register int i,j,t;
register site *st;
msg_tag *tag;
complex sum;
complex plp;

    sum = cmplx(0.0,0.0);
    FORALLSITES(i,st){lattice[i].tempmat2 = lattice[i].link[TUP];}
    for(t=1;t<nt;t++){
	tag=start_gather( F_OFFSET(tempmat2), sizeof(su3_matrix),
	    TUP, EVENANDODD, gen_pt[0] );
	wait_gather(tag);
        FORALLSITES(i,st){
	    if( st->t != 0 )continue;
	    if(t==1){
		mult_su3_nn( &(st->link[TUP]), gen_pt[0][i],
		    &(st->staple));
	    }
	    else {
		mult_su3_nn( &(st->staple), gen_pt[0][i],
		    &(st->tempmat2));
		lattice[i].staple = lattice[i].tempmat2;
	    }
	}
        FORALLSITES(i,st){
	    st->tempmat1 = *(su3_matrix *)(gen_pt[0][i]);
	}
        FORALLSITES(i,st){
	    st->tempmat2 = st->tempmat1;
	}
	cleanup_gather(tag);
    }
    FORALLSITES(i,st){
	if( st->t != 0 )continue;
	plp = trace_su3( &(st->staple) );
	CSUM(sum,plp);
    }
    g_complexsum( &sum );
    plp.real = sum.real /((float)(nx*ny*nz));
    plp.imag = sum.imag /((float)(nx*ny*nz));
    return(plp);
}
