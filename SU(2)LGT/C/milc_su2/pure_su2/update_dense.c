/********** update_dense.c ****************************************************/
/* MIMD version 4 */

/*
 Update lattice with Metropolis algorithm, with "almost quenched"
 high density.
*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

int update_dense()  {
   int Nhit;
   int step, iters=0;
   void relax_space(),monte_space(),monte_time(),reunitarize(),check_unitarity();

    /* check unitarity before doing anything */
	check_unitarity();
	
    /* do "steps" overrelaxed steps and stepsQ  qhb steps */
	relax_space(steps); 
	monte_space(stepsQ);
	monte_time(stepsQ);

        /* reunitarize the gauge field */

        reunitarize();         

    if(steps > 0)return (iters/steps);
    else return(-99);
}
