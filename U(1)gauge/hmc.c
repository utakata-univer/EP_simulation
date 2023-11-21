#include "u1.h"

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
static double mom[L+2][L+2][L+2][L+2][D];
/* saved phi field for hmc */
static double phiold[V];
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
    int nx,ny,nz,nt,mu;

    /* H=p^2/+S*/

    act=action();

    for (nx=0;nx<L;nx++)
    {
    for (ny=0;ny<L;ny++)
    {
    for (nz=0;nz<L;nz++)
    {
    for (nt=0;nt<L;nt++)
    {
    	 for (mu=0;mu<D;mu++)
	 {      
		 act+=mom[nx][ny][nz][nt][mu]*mom[nx][ny][nz][nt][mu]/2;
         }
    }
    }
    }
    }
    return act;
}

{
    a[0][ny][nz][nt][mu]=a[L+1][ny][nz][nt][mu];
    a[1][ny][nz][nt][mu]=a[L][ny][nz][nt][mu];
}

void initialize();
{
     
     
     ranlxd(b,)

     a[]=b

}
/********************************************************************** 
 *     move_phi()
 *     elementary leap frog step for the update of the phi field
 *     does phi <- phi+mom*eps
 **********************************************************************/
static void move_phi(double eps)
{
    int nx,ny,nz,nt,mu;

    for (nx=0;nx<L;nx++)
    {
    for (ny=0;ny<L;ny++)
    {
    for (nz=0;nz<L;nz++)
    {
    for (nt=0;nt<L;nt++)
    {
         for (mu=0;mu<D;mu++)
         {
              a[nx][ny][nz][nt][mu]+=mom[nx][ny][nz][nt][mu]*eps;
         }
    }
    }
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
    int nx,ny,nz,nt,mu,nu;
    double  force;


    for (i=0;i<V;i++)
    { 
       force=0; 
    
       for (nu=0;nu<4;nu++);

       force+=2*sin(a[nx][ny][nz][nt][mu]+a[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]-a[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-a[nx][ny][nz][nt][nu])+2*sin(a[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]+a[nx-d[0][nu]+d[0][mu]][ny-d[1][nu]+d[1][mu]][nz-d[2][nu]+d[2][mu]][nt-d[2][nu]+d[3][mu]][nu]-a[nx][ny][nz][nt][mu]-a[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]);
	mom[nx][ny][nz][nt][mu]+=force*eps;
    }

}

/********************************************************************** 
 * update_hmc()
 * Do one trajectory with of the HMC
 **********************************************************************/

static void update_hmc(hmc_params_t* pars)
{
    int i, istep, nstep=pars->nstep;
    double eps=(pars->tlength)/(pars->nstep);
    double startH, endH, deltaH;
    double r;
    double mom2[L*L*L*L*L+L*L*L*L+L*L*L+L*L+4];

    /*
     *  Action: Sum_x -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2
     */

    gauss_rand(L*L*L*L*L+L*L*L*L+L*L*L+L*L+4,mom2);

    /* refresh the momenta */
    for (nx=1;nx<L+1;nx++)
   {
    for (ny=1;ny<L+1;ny++)
   {
    for (nz=1;nz<L+1;nz++)
   {
    for (nt=1;nt<L+1;nt++)
   {	    	  
    for (mu=1;mu<L+1;mu++)
   {	    
	   mom2[L*L*L*L*nx+L*L*L*ny+L*L*nz+L*nt+mu]=mom[nx][ny][nz][nt][mu];
   }

   }

   }

   }

   } 
    /* keep old phi field */
    for (nx=1;nx<L+1;nx++)
   {
    for (ny=1;ny<L+1;ny++)
   {
    for (nz=1;nz<L+1;nz++)
   {
    for (nt=1;nt<L+1;nt++)
   {
    for (mu=0;mu<4;mu++)
   {
     aold[nx][ny][nz][nt][mu]=a[nx][ny][nz][nt][mu];
   }

   }

   }

   }

   } 
    /* measure hamiltonian */
    startH=hamilton();

    /* do the trajectory */
    for (istep=0;istep<nstep;istep++)
    {
	move_a(eps/2.);
	move_m  (eps);
	move_a(eps/2.);
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
        	    for (nx=1;nx<L+1;nx++)
                 {                    for (ny=1;ny<L+1;ny++)
                 {
                    for (nz=1;nz<L+1;nz++)
                 {
                    for (nt=1;nt<L+1;nt++)
                 {
                    for (mu=1;mu<L+1;mu++)
                 {
	            a[nx][ny][nz][nt][mu]=aold[nx][ny][nz][nt][mu];
                 }
                 }
		 }
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

    g=apars->g

    accept=reject=0;

    expdH=0.0;
    for (isweep=0;isweep<hpars->ntraj;isweep++)
    {
	update_hmc(hpars);
	measure();
	if((isweep+1)%(hpars->naccu)==0) {
	    print_meas();
	    printf("EXPDH %e\n",expdH/hpars->naccu);
	    expdH=0.0;
	}
    }

    return ((double)accept)/((double)accept+reject);
}
