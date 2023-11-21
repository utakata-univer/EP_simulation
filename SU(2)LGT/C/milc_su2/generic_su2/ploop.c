/****************** ploop.c ************************************/
/* MIMD version 3 */
/* evaluate the Polyakov loops */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include LATDEF
#include <comdefs.h>

complex ploop() {
register int i,j;
register site *st;
su3_matrix mat1,mat2;
register su3_matrix *mp1,*mp2,*mpt;
complex sum;
complex plp;
int x,y,z,t;

    sum = cmplx(0.0,0.0);
    FORALLSITES(i,st){
	if( st->t != 0 )continue;
	x=st->x; y=st->y; z=st->z; t=st->t;
	mp1 = &(st->link[TUP]);
	mp2 = &mat1;
	for (j=1;j< nt;j++) {
	    neighbor_coords(x,y,z,t,TUP,&x,&y,&z,&t);
	    mpt = (su3_matrix *)field_pointer_at_coordinates(
		F_OFFSET(link[TUP]), sizeof(su3_matrix), x,y,z,t);
	    mult_su3_nn( mp1 , mpt , mp2 );
	    cleanup_field_pointer(mpt);
	    mpt = mp1;
	    mp1 = mp2;
	    if(j==1)mp2 = &mat2; else mp2 = mpt;
	}
	plp = trace_su3(mp1);
	CSUM(sum,plp);
    }
    g_complexsum( &sum );
    plp.real = sum.real /((float)(nx*ny*nz));
    plp.imag = sum.imag /((float)(nx*ny*nz));
    return(plp);
}
