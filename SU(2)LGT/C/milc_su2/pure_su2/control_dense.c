/**************************** control_dense.c ******************************/
/* Main procedure for pure gauge SU2 */
/* Almost quenched at nonzero density */
/* MIMD version 4 */

/* We use over-relaxed + heat bath updatings in this code. */

#define CONTROL
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"    /* global variables for lattice fields */
#include <comdefs.h>

main(argc,argv)  {
   int readin();
   int meascount,todo;
   int prompt;
   float ssplaq,stplaq,rpbpp,rpbpm;
   float rsq;
   complex plp,ploop();
   void plaquette(),save_ascii();
   double dtime;
   int setup();
   int i;

   initialize_machine(argc,argv);

   g_sync();
   /* set up */
   prompt = setup();

   /* loop over input sets */
   while( readin(prompt) == 0){

        /* perform warmup trajectories */
        dtime = -dclock();
 
 /* call plaquette measuring process */
                plaquette(&ssplaq,&stplaq);
                if(this_node==0)printf("START %e %e\n",
                    (double)ssplaq,(double)stplaq);

        for(todo=warms; todo > 0; --todo ){
            update_dense();
        }
        if(this_node==0)printf("THERMALIZATION COMPLETED\n");

        /* perform measuring trajectories, reunitarizing and measuring  */
        meascount=0;            /* number of measurements               */
        plp = cmplx(99.9,99.9);
        for(todo=trajecs; todo > 0; --todo ){ 

            /* do the trajectories */
            update_dense();

            /* measure every "propinterval" trajectories */
            if((todo%propinterval) == 0){
            
                /* call plaquette measuring process */
                plaquette(&ssplaq,&stplaq);

                /* don't bother to */
                /* call the Polyakov loop measuring program */
                /* plp = ploop(); */

                ++meascount;
                if(this_node==0)printf("GMES %e %e %e %e %e\n",
                    (double)plp.real,(double)plp.imag,99.9,
                    (double)ssplaq,(double)stplaq);
                /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */

                fflush(stdout);
            }
        }       /* end loop over trajectories */

        if(this_node==0)printf("RUNNING COMPLETED\n");

        dtime += dclock();
        if(this_node==0){
            printf("Time = %e seconds\n",dtime);
        }
        fflush(stdout);
	dtime = -dclock();

        /* save lattice if requested */
        if( saveflag == SAVE_ASCII){
            save_ascii(savefile,beta,0.0);
            dtime += dclock();
            if(this_node==0){
                printf("Time for saving lattice = %e seconds\n",dtime);
            }
            fflush(stdout);
        }
        if( saveflag == SAVE_BINARY){
            save_binary(savefile,beta,0.0);
            dtime += dclock();
            if(this_node==0){
                printf("Time for saving lattice = %e seconds\n",dtime);
            }
            fflush(stdout);
        }
    }
}

/* --------- converted to su2  --------*/




