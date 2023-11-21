#include "su2.h"

/*
 *  File heatbath.c
 *  Heat Bath algorithm for the su2 gauge theory
 *
 *  double heatbath(act_params_t *apars, heatbath_params_t *hpars)
 *
 *  static void update_heatbath(hmc_params_t* pars)
 *
 */

static double beta;

/********************************************************************** 
   update_heatbath()
   Do one sweep of the Heat Bath Update 

   -S           = -Sum_{x, mu<nu} beta* ReTr(1-Up)/Nc

                (For SU(2), Nc=2 and Tr(U) is real)

   -dS[U(x,mu)] = beta* Tr( U(x,mu) Stpl(x,mu) ) 
                
   Stpl(x,mu) =sum_{nu \= mu} [ Stpl+(x,nu) + Stpl_(x,nu) ]
            
  Stpl+(x,mu)=U(x+\hat nu) {U(x,nu) U(x+\hat nu, mu)}^-1
  Stpl-(x,mu)={U(x-\hat nu,mu) U(x-\hat nu+\hat mu, nu)}^-1 U(x-\hat \nu,nu)

  Stpl = k \bar{U}

  r1, r2, r3, r4: uniform random number \in [0,1];

  a0 = 1+ (1/beta k) log ( exp(-2 \beta k) + (1-exp(-2 \beta k)) r1 ) 

  accept if 0 < r2 < \sqrt(1-a0^2) 

  c = -1 + 2*r3
  s = sqrt(1-c*c)

  a1 = \sqrt(1-a0^2) s cos (2\pi r4) 
  a2 = \sqrt(1-a0^2) s sin (2\pi r4) 
  a3 = \sqrt(1-a0^2) c
    
  update U  to  (a0 + i \sigma_k a_k) \bar{U}^{-1} 

 **********************************************************************/

static void update_heatbath()
{
  int n,m,l,mu,nu;  
  double s0,s1,s2,s3;
  double k,bk;
  double p0,p1,p2,p3;
  double m0,m1,m2,m3;

  double u0,u1,u2,u3;
  double v0,v1,v2,v3;

  double r1,r2,r3,r4;

  double Pi;
  double r,c,s,phi;

  int ihit;
  /* double P; */

  Pi=4.* atan(1.);

  for(n=0;n<V;n++){
    for(mu=0;mu<D;mu++){
	
      /* compute the sum of the staples */
      s0=0.0;
      s1=0.0;
      s2=0.0;
      s3=0.0;

      for(nu=0;nu<D;nu++){
	if(nu != mu){
	  /* printf(" %i %i \n", mu,nu); */

	  /* upward staple */
	  m=hop[n][mu]; 
	  u0= a0[m][nu];
	  u1= a1[m][nu];
	  u2= a2[m][nu];
	  u3= a3[m][nu];

	  m=hop[n][nu]; 
	  v0= a0[n][nu]*a0[m][mu] - a1[n][nu]*a1[m][mu] - a2[n][nu]*a2[m][mu] - a3[n][nu]*a3[m][mu];
	  v1= a0[n][nu]*a1[m][mu] + a1[n][nu]*a0[m][mu] - a2[n][nu]*a3[m][mu] + a3[n][nu]*a2[m][mu];
	  v2= a0[n][nu]*a2[m][mu] + a2[n][nu]*a0[m][mu] - a3[n][nu]*a1[m][mu] + a1[n][nu]*a3[m][mu];
	  v3= a0[n][nu]*a3[m][mu] + a3[n][nu]*a0[m][mu] - a1[n][nu]*a2[m][mu] + a2[n][nu]*a1[m][mu];

	  p0= u0*v0+u1*v1+u2*v2+u3*v3;
	  p1=-u0*v1+u1*v0+u2*v3-u3*v2;
	  p2=-u0*v2+u2*v0+u3*v1-u1*v3;
	  p3=-u0*v3+u3*v0+u1*v2-u2*v1;

	  /* printf("up staple for %i %i %i : %e %e %e %e\n",n,mu,nu,p0,p1,p2,p3); */

	  s0+=p0;
	  s1+=p1;
	  s2+=p2;
	  s3+=p3;

	  /* downward staple */
	  l=hop[n][nu+D];

	  m=hop[l][mu]; 
	  u0= a0[l][mu]*a0[m][nu] - a1[l][mu]*a1[m][nu] - a2[l][mu]*a2[m][nu] - a3[l][mu]*a3[m][nu];

	  u1= a0[l][mu]*a1[m][nu] + a1[l][mu]*a0[m][nu] - a2[l][mu]*a3[m][nu] + a3[l][mu]*a2[m][nu];
	  u2= a0[l][mu]*a2[m][nu] + a2[l][mu]*a0[m][nu] - a3[l][mu]*a1[m][nu] + a1[l][mu]*a3[m][nu];
	  u3= a0[l][mu]*a3[m][nu] + a3[l][mu]*a0[m][nu] - a1[l][mu]*a2[m][nu] + a2[l][mu]*a1[m][nu];

	  v0= a0[l][nu];
	  v1= a1[l][nu];
	  v2= a2[l][nu];
	  v3= a3[l][nu];

	  p0= u0*v0+u1*v1+u2*v2+u3*v3;
	  p1= u0*v1-u1*v0+u2*v3-u3*v2;
	  p2= u0*v2-u2*v0+u3*v1-u1*v3;
	  p3= u0*v3-u3*v0+u1*v2-u2*v1;

	  /* printf("dw staple for %i %i %i : %e %e %e %e\n",n,mu,nu,p0,p1,p2,p3); */

	  s0+=p0;
	  s1+=p1;
	  s2+=p2;
	  s3+=p3;

	}
      }

      /* printf("staple for %i %i    : %e %e %e %e\n",n,mu,s0,s1,s2,s3); */

      /* compute the factor k */
      k=sqrt(s0*s0+s1*s1+s2*s2+s3*s3);

      /* compute \bar U^{-1} */
      m0 = s0/k;
      m1 =-s1/k;
      m2 =-s2/k;
      m3 =-s3/k;

      bk = beta*k;

      /* printf("staple for %i %i    : %e %e %e %e %e %e %e\n",n,mu,m0,-m1,-m2,-m3,k,bk,beta); */

      /* generate U' = U \bar U */
      ihit=0;
      do{
	ranlxd(&r1,1);
	p0 = 1.0 + log( exp(-2.0*bk)+(1.0-exp(-2.0*bk))*r1 )/bk;
	r = sqrt(1-p0*p0);
	ranlxd(&r2,1);
        /* printf("%e %e %e\n",r1,r2,r); */
	ihit++;
      }while(r < r2);
      /* printf("ihit %i\n",ihit); */

      ranlxd(&r3,1);
      c = -1.0 +2.0*r3;
      s = sqrt(1.0-c*c);

      ranlxd(&r4,1);
      phi =2.0*Pi*r4;
        
      p1=r*s*cos(phi);
      p2=r*s*sin(phi);
      p3=r*c;

      /* printf("U' for %i %i    : %e %e %e %e\n",n,mu,p0,p1,p2,p3); */

      /* compute U = U' \bar U^{-1} */
      a0[n][mu]=p0*m0-p1*m1-p2*m2-p3*m3;
      a1[n][mu]=p1*m0+p0*m1-p2*m3+p3*m2;
      a2[n][mu]=p2*m0+p0*m2-p3*m1+p1*m3;
      a3[n][mu]=p3*m0+p0*m3-p1*m2+p2*m1;

      /*
      for(m=0;m<V;m++){
	for(nu=0;nu<D;nu++){
	  printf("U(%i,%i)= %e %e %e %e\n",m,nu,a0[m][nu],a1[m][nu],a2[m][nu],a3[m][nu]);
	}
      }
      P=plaquette();
      */

    }
  }

}


/********************************************************************** 
 * heatbath()
 * Does n sweeps of the Heat Bath MC algorithm. Measurement after each
 * sweep, the averaged measured values are printed out in fixed
 * intervals.
 **********************************************************************/

double heatbath(act_params_t *apars, heatbath_params_t *hpars)
{
  int isweep, itherm;
  double P;

  beta=apars->beta;


  for (itherm=0;itherm<hpars->ntherm;itherm++)
    {
      update_heatbath();

      P=plaquette();
      printf("%i %e\n",itherm+1,1.0-P);

    }

/*  init_measure();*/

  for (isweep=0;isweep<hpars->nsweep;isweep++)
    {
      update_heatbath();

/*      P=plaquette();
      printf("%i %e\n",isweep,P); */
     

	if((isweep+1)%(hpars->naccu)==0){
/*          measure();*/
	  P=plaquette();
          printf("%e\n",-log(P));	
	}
                
    }
/*     print_meas(); */

  return 0;
}
