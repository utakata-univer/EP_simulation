/****************** ploop3.c ************************************/
/* MIMD version 3 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */

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
int d[4];

    if(ZSKEW != 2){printf("PLOOP: FIX ME!!\n"); terminate(0);}

    sum = cmplx(0.0,0.0);
    d[XUP] = d[YUP] = d[ZUP] = 0;
    /* First multiply the link on every even site by the link above it */
    /* We will compute the Polyakov loop "at" the even sites in the 
	first two time slices. */
    tag=start_gather( F_OFFSET(link[TUP]), sizeof(su3_matrix),
	TUP, EVEN, gen_pt[0] );
    wait_gather(tag);
    FOREVENSITES(i,st){
	mult_su3_nn( &(st->link[TUP]), gen_pt[0][i], &(st->tempmat1));
    }
    cleanup_gather(tag);

    for(t=2;t<nt;t+=2){
	d[TUP] = t;	/* distance from which to gather */
	tag=start_general_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
	    d, EVEN, gen_pt[0] );
	wait_general_gather(tag);
        FOREVENSITES(i,st){
	    if( st->t > 1 )continue;  /* only compute on first two slices */
	    mult_su3_nn( &(st->tempmat1), gen_pt[0][i], &(st->tempmat2));
	    lattice[i].tempmat1 = lattice[i].tempmat2;
	    /* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
	}
	cleanup_general_gather(tag);
    }

    tag=start_gather( F_OFFSET(link[ZUP]), sizeof(su3_matrix),
	ZUP, EVEN, gen_pt[0] );
    wait_gather(tag);
    FOREVENSITES(i,st){
	if( st->t > 1 )continue;  /* only compute on first two slices */
	mult_su3_nn( &(st->link[ZUP]), gen_pt[0][i], &(st->tempmat2));
	mult_su3_na( &(st->tempmat1), &(st->tempmat2), &(st->staple) );
	lattice[i].tempmat1 = lattice[i].staple;
    }
    cleanup_gather(tag);

    FOREVENSITES(i,st){
	if( st->t > 1 )continue;
	plp = trace_su3( &(st->tempmat1) );
	CSUM(sum,plp);
    }
    g_complexsum( &sum );
    plp.real = sum.real /((float)(nx*ny*nz));
    plp.imag = sum.imag /((float)(nx*ny*nz));
    return(plp);
}
