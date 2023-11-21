/****** dsdu_qhb.c  -- compute the staple ******************/

/* MIMD version 3 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

void
dsdu_qhb(dir1,parity) 
register int dir1;
int parity;
{
register int i,dir2, otherparity;
register site *st;
msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
int start, first,count;
su2_matrix tmat1,tmat2;
int disp[4];	/* displacement vector for general gather */
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	start=1; /* indicates staple sum not initialized */

	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1)
	{
	    /* displacement vector for link 2 sites away */
	    for(i=XUP;i<=TUP;i++)disp[i]=0;
	    disp[dir1] = 1;
	    disp[dir2] = -1;

	    /* get link[dir2] from direction dir1 */
	    tag0 = start_gather( F_OFFSET(link[dir2]), sizeof(su2_matrix),
		dir1, parity, gen_pt[0] );

	    /* get link[dir1] from direction dir2 */
	    tag1 = start_gather( F_OFFSET(link[dir1]), sizeof(su2_matrix),
		dir2, parity, gen_pt[1] );

	    /* get link[dir2] from direction -dir2 */
	    tag2 = start_gather( F_OFFSET(link[dir2]), sizeof(su2_matrix),
		OPP_DIR(dir2), parity, gen_pt[2] );

	    /* get link[dir1] from direction -dir2 */
	    tag3 = start_gather( F_OFFSET(link[dir1]), sizeof(su2_matrix),
		OPP_DIR(dir2), parity, gen_pt[3] );

	    /* get link[dir2] from displacement +dir1-dir2 */
	    tag4 = start_general_gather( F_OFFSET(link[dir2]),
		sizeof(su2_matrix), disp, parity, gen_pt[4] );

	    /* Upper staple */
	    wait_gather(tag0);
	    wait_gather(tag1);
          if(start){  /* this is the first contribution to staple */
	FORSOMEPARITY(i,st,parity){
	        mult_su2_nn( &(st->link[dir2]), 
			    (su2_matrix *)gen_pt[1][i], &tmat1 );
		mult_su2_na( &tmat1, (su2_matrix *)gen_pt[0][i], &(st->staple) );
		}
		start=0; 
		}
		else{
	FORSOMEPARITY(i,st,parity){
		mult_su2_nn( &(st->link[dir2]), 
			    (su2_matrix *)gen_pt[1][i], &tmat1 );
		mult_su2_na( &tmat1, (su2_matrix *)gen_pt[0][i], &tmat2 );
		add_su2_matrix( &(st->staple), &tmat2, &(st->staple));
		}
	    } /* upper staple */
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);

	    /* Lower staple */
	    wait_gather(tag2);
	    wait_gather(tag3);
	    wait_general_gather(tag4);
	FORSOMEPARITY(i,st,parity){
	        mult_su2_an( (su2_matrix *)gen_pt[2][i], 
			    (su2_matrix *)gen_pt[3][i], &tmat1 );
	        mult_su2_nn( &tmat1, (su2_matrix *)gen_pt[4][i], &tmat2 );
		add_su2_matrix( &(st->staple), &tmat2, &(st->staple));
	    }  /* lower staple */
	    cleanup_gather(tag2);
	    cleanup_gather(tag3);
	    cleanup_general_gather(tag4);
	}
}

/*su2 */











