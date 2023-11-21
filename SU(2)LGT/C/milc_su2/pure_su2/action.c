/*************** action.c ****************************************/
/* MIMD version 3 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

float action(){
void plaquette();
float hmom_action();
float ssplaq,stplaq,g_action,h_action;
    plaquette(&ssplaq,&stplaq);
    g_action = -beta*volume*(ssplaq+stplaq);
    h_action = hmom_action();
if(this_node==0)printf("ACTION: g,h = %e  %e  %e\n",
(double)g_action,(double)h_action, (double)(g_action+h_action));
    return(g_action+h_action);
}

/* gauge momentum contribution to the action */
float hmom_action() {
register int i,dir;
register site *s;
float sum;
float ahmat_mag_sq();

    sum=0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
            sum += ahmat_mag_sq( &(s->mom[dir]) );
	}
    }
    g_floatsum( &sum );
    return(sum);
}

/* magnitude squared of an antihermition matrix */
float ahmat_mag_sq(pt) anti_hermitmat *pt; {
register float x,sum;
    x = pt->m00im; sum  = 0.5*x*x;
    x = pt->m11im; sum += 0.5*x*x;
    x = pt->m22im; sum += 0.5*x*x;
    x = pt->m01.real; sum += x*x;
    x = pt->m01.imag; sum += x*x;
    x = pt->m02.real; sum += x*x;
    x = pt->m02.imag; sum += x*x;
    x = pt->m12.real; sum += x*x;
    x = pt->m12.imag; sum += x*x;
    return(sum);
}

/* copy a gauge field - an array of four su2_matrices */
gauge_field_copy(src,dest) field_offset src,dest; {
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    su2mat_copy( F_PT(s,src2), F_PT(s,dest2) );
	    src2 += sizeof(su2_matrix);
	    dest2 += sizeof(su2_matrix);
	}
    }
}
