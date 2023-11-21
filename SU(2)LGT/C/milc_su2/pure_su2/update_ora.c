/********** update_ora.c ****************************************************/
/* MIMD version 4 */

/*
 Update lattice with microcanonical overrelaxed algorithm
*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

int update()  {
   int Nhit;
   int step, iters=0;
   void relax(),monte(),reunitarize(),check_unitarity();

    /* check unitarity before doing anything */
	check_unitarity();
	
    /* do "steps" overrelaxed steps and stepsQ  qhb steps */
	relax(steps); 
	monte(stepsQ);

        /* reunitarize the gauge field */

        reunitarize();         

    if(steps > 0)return (iters/steps);
    else return(-99);
}
