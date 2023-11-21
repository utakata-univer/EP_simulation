/****************** ploop3.c ************************************/
/* MIMD version 3 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */

#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"
#include LATDEF
#include <comdefs.h>

float ploop() {
register int i,j,t;
register site *st;
msg_tag *tag;
float sum, plp;
int d[4];

    sum = 0.0;
    d[XUP] = d[YUP] = d[ZUP] = 0;
    /* First multiply the link on every even site by the link above it */
    /* We will compute the Polyakov loop "at" the even sites in the 
	first two time slices. */
    tag=start_gather( F_OFFSET(link[TUP]), sizeof(su2_matrix),
	TUP, EVEN, gen_pt[0] );
    wait_gather(tag);
    FOREVENSITES(i,st){
	mult_su2_nn( &(st->link[TUP]), (su2_matrix *)gen_pt[0][i], 
		    &(st->tempmat1));
    }
    cleanup_gather(tag);

    for(t=2;t<nt;t+=2){
	d[TUP] = t;	/* distance from which to gather */
	tag=start_general_gather( F_OFFSET(tempmat1), sizeof(su2_matrix),
	    d, EVEN, gen_pt[0] );
	wait_general_gather(tag);
        FOREVENSITES(i,st){
	    if( st->t > 1 )continue;  /* only compute on first two slices */
	    mult_su2_nn( &(st->tempmat1), (su2_matrix *)gen_pt[0][i], 
			&(st->tempmat2));
	    lattice[i].tempmat1 = lattice[i].tempmat2;
	    /* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
	}
	cleanup_general_gather(tag);
    }
    FOREVENSITES(i,st){
	if( st->t > 1 )continue;
	plp = trace_su2( &(st->tempmat1) );
	/* Save for later correlation measurements */
	sum += plp;
#ifdef BPCORR
	/* Save for subsequent correlation measurements */
	/* Note the results are saved on even sites in
	   slices 0 and 1 */
	st->ploop = plp;
#endif
    }
    g_floatsum( &sum );
    plp = sum /((float)(nx*ny*nz));
    return(plp);
}
