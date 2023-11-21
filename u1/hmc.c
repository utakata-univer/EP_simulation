#include "phi4.h"

/*
 *  File hmc.c
 *  Molecular dynamics algorithm for the phi4 theory
 *
 *  double hmc(act_params_t *apars, hmc_params_t *hpars)
 *       Hybrid Monte Carlo algorithm, starts from the global field phi
 *       The parameters of the action (kappa and lambda) are passed by
 *       the apars, the parameters of the algorithm by hpars, here
 *       we need the number of trajectories (ntraj), the trajectory length
 *       (tlength) and the number of steps per trajectory (nstep).
 *       calls update_md() ntraj times.
 *       Returns the acceptance rate.
 *
 *
 *  static double hamilton()
 *       computes the value of the HMC Hamiltonian H=mom^2/2+S(phi)
 *
 *  static void move_phi(double eps)
 *       one of the two elementary building blocks of the leapfrog
 *       phi <- phi + eps*mom
 *
 *  static void move_m(double eps)
 *       the other elementary building block of the leapfrog
 *       mom <- mom - eps*dS/dphi(phi)
 *
 *
 *  static void update_hmc(hmc_params_t* pars)
 *       The actual workhorse for the hmc routine.
 *       Does one trajectory: 
 *       momentum heatbath, leapfrog integration and 
 *       acceptance step
 *
 */



/* HMC momenta */
static double mom[V][4];
/* saved phi field for hmc */
static double aold[V][4];
/* book keeping for acceptance rate */
static int accept, reject;
static double expdH;

static double g;

/**********************************************************************
 *     hamilton
 **********************************************************************/


static double hamilton()
{
    double act;
    int i,mu;

    /* H=p^2/+S*/

    act=action();

    for (i=0;i<V;i++) {
	   for (mu=0;mu<D;mu++) {
		   act+=mom[i][mu]*mom[i][mu]/2;
	   }
                      }
    return act;
}


/********************************************************************** 
 *     move_phi()
 *     elementary leap frog step for the update of the phi field
 *     does phi <- phi+mom*eps
 **********************************************************************/
static void move_phi(double eps)
{
    int i,mu;

    for (i=0;i<V;i++)
    {
	for (mu=0;mu<D;mu++) {
		a[i][mu]+=mom[i][mu]*eps;
	}
    }

}

/********************************************************************** 
 *   move_m()
 *   elementary leap frog step for the update of the momenta
 *   does mom <- mom-eps*force
 **********************************************************************/

static void move_m( double eps)
{
    int i,mu,nu;
    double force;


    for (i=0;i<V;i++)
    {
	for (mu=0;mu<D;mu++)
	{   
	force=0;       
		for (nu=0;nu<D;nu++){
                                         
					 if (nu != mu){

	force+=(2/g*g)*(sin(a[i][mu]+a[hop[i][mu]][nu]-a[hop[i][nu]][mu]-a[i][nu])+sin(a[i][mu]-a[hop[hop[i][mu]][D+nu]][nu]-a[hop[i][D+nu]][mu]+a[hop[i][D+nu]][nu]));

					 }
                                
	                        	}
		 mom[i][mu]+=force*eps;
				 
         }
    }
	
}

/********************************************************************** 
 * update_hmc()
 * Do one trajectory with of the HMC
 **********************************************************************/

static void update_hmc(hmc_params_t* pars)
{
    int i,mu, istep, nstep=pars->nstep;
    double eps=(pars->tlength)/(pars->nstep);
    double startH, endH, deltaH;
    double r,mom2[D*V];

    /*
     *  Action: Sum_x -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2
     */

    /* refresh the momenta */
    gauss_rand(D*V,mom2);

    for (i=0;i<V;i++){
           for (mu=0;mu<D;mu++) {
	    mom[i][mu]=mom2[i*D+mu];
	   }
	             }

    /* keep old phi field */
    for (i=0;i<V;i++){
	for (mu=0;mu<D;mu++) {
	     	aold[i][mu]=a[i][mu];
	}
	             }

    /* measure hamiltonian */
    startH=hamilton();

    /* do the trajectory */
    for (istep=0;istep<nstep;istep++)
    {
	move_phi(eps/2.);
	move_m  (eps);
	move_phi(eps/2.);
    }

    /* compute energy violation */
    endH=hamilton();
    deltaH=endH-startH;
    expdH+=exp(-deltaH); 

    /* acceptance step */
    if (deltaH<0) {
	accept++;
    } else {
	ranlxd(&r,1);
	if (exp(-deltaH)>r) {
	    accept++;
	} else {
	    reject++;
	    for (i=0;i<V;i++){
		for (mu=0;mu<D;mu++) {
		     	a[i][mu]=aold[i][mu];
		}
	}
    }



}
}

/********************************************************************** 
 * hmc()
 * Does n trajectories of the HMC algorithm. Measurement after each
 * trajectory, the averaged measured values are printed out in fixed
 * intervals.
 **********************************************************************/

double hmc(act_params_t *apars, hmc_params_t *hpars)
{
    int isweep;

    g=apars->g;
	    
    accept=reject=0;

    expdH=0.0;
    for (isweep=0;isweep<hpars->ntraj;isweep++)
    {
	update_hmc(hpars);
/*	measure();*/
	if((isweep+1)%(hpars->naccu)==0) {
/*	    print_meas();*/
	    printf("EXPDH %e\n",expdH/hpars->naccu);
	    expdH=0.0;
	}
    }

    return ((double)accept)/((double)accept+reject);
}

