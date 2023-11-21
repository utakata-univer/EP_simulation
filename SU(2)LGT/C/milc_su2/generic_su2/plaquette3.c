/************************** plaquette3.c *******************************/
/* MIMD version 3 */
/* This version uses gathers to get the neighbors */
/* version of 8/18/91 by CWB to use less memory: 
    it's ok to store things in staple since that is needed
    anyway by update, but all use of fields tempmat1 and tempmat2
    is eliminated    */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include LATDEF /* global variables for lattice fields */
#include <comdefs.h>

void plaquette(ss_plaq,st_plaq) float *ss_plaq,*st_plaq; {
register int i,dir1,dir2;
register site *s;
register su2_matrix *m1,*m4;
su2_matrix mtmp;
float ss_sum,st_sum;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather( F_OFFSET(link[dir2]), sizeof(su2_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather( F_OFFSET(link[dir1]), sizeof(su2_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->link[dir1]);
		m4 = &(s->link[dir2]);
		mult_su2_an(m4,m1,&(s->staple));
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su2_nn( &(s->staple),(su2_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP )st_sum +=
		    realtrace_su2((su2_matrix *)(gen_pt[1][i]),&mtmp);
		else          ss_sum +=
		    realtrace_su2((su2_matrix *)(gen_pt[1][i]),&mtmp);
	    }

	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    g_floatsum( &ss_sum );
    g_floatsum( &st_sum );
    *ss_plaq = ss_sum /((float)(3*nx*ny*nz*nt));
    *st_plaq = st_sum /((float)(3*nx*ny*nz*nt));
} /* plaquette3 */

/* converted to su2 */
