/************************** monte.c *******************************/
/* Kennedy-Pendleton quasi heat bath on SU(2) subgroups */
/* MIMD version 3 */
/* T. DeGrand March 1991 */
/* J. Hetrick and D. Toussaint April 1995 */


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

/* Generic definitions - could be useful elsewhere */
#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/


void monte(int NumStp) {
/* Do K-P quasi-heat bath*/
  int NumTrj,Nhit;
  int parity;
  float xr1,xr2,xr3,xr4;
  float r,r2,rho,z;
  float al,d, xl,xd;
  float a0,a1,a2,a3;
  int  k,kp, cr, nacd, test;
  float avekp, avecr;
  float pi2, b3;
  int first,count;
  register int dir,i;
  register site *st;
  su2_matrix action, h, a;
  pi2= 2.0*PI;
  b3=beta/3.0;

/* fix bug by adding loop over NumTrj; before 1 (and only 1) heat bath
   hit was dome, regardless of NumStp    */
  for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++) {

/* fix bug by looping over odd AND even parity */
    for(parity=ODD;parity<=EVEN;parity++) {
      FORALLUPDIR(dir) {
/* compute the gauge force */
	dsdu_qhb(dir,parity); 
/* now for the qhb updating */
	FORSOMEPARITY(i,st,parity) {
	  z = sqrt((double)det_su2(&(st->staple)));
	  reunit_su2(&(st->staple));
	  mult_su2_na(&(st->link[dir]), &(st->staple), &action);

printf("start site %d, z = %f\n",i,z);
dump_su2_mat(&(st->link[dir]));
dump_su2_mat(&(st->staple));

/* generate gaussian random number about the Identity of SU(2)

/* get four random numbers (add a small increment to prevent taking log(0.)*/
	  xr1=myrand(&(st->site_prn));
	  xr1 = (log((double)(xr1+ 1.e-10)));
	  xr2=myrand(&(st->site_prn));
	  xr2 = (log((double)(xr2+ 1.e-10)));
	  xr3=myrand(&(st->site_prn));
	  xr4=myrand(&(st->site_prn));

/* debug:
if(this_node == 0)printf("rand= %e %e %e %e\n",xr1,xr2,xr3,xr4); 
*/
/* generate an su(2) matrix h according to exp(bg/3 * re tr(h*s))
   rewrite re tr(h*s) as re tr(h*v)z where v is
   an su(2) matrix and z is a real normalization constant.
   Let v = z*v. (z is 2*xi in k-p notation)
   v is represented in the form v(0) + i*sig*v (sig are pauli)
   v(0) and vector v are real.
   Let a = h*v and now generate a
   rewrite beta/3 * re tr(h*v) * z as al*a0
   a0 has prob(a0) = n0 * sqrt(1 - a0**2) * exp(al * a0)
*/
	  al=b3*z;

/* debug:
if(this_node == 0)printf("al= %e\n",al);
*/

/* let a0 = 1 - del**2, get d = del**2
   such that prob2(del) = n1 * del**2 * exp(-al*del**2) */

	  d= -(xr2  + xr1*xr3*xr3)/al;

/* monte carlo prob1(del) = n2 * sqrt(1 - 0.5*del**2)
   then prob(a0) = n3 * prob1(a0)*prob2(a0) */

/* now  beat each  site into submission */
	  nacd = 0;
	  if ((1.00 - 0.5*d) > xr4*xr4) nacd=1;

/*------ k-p algorithm --------*/	  
	  if(nacd == 0 && al > 2.0) {  
printf("KP method\n");
	    test=0;
	    for(k=0;k<20 && test == 0;k++) {
	      kp++;
/*  get four random numbers (add a small increment to prevent taking log(0.)*/
	      xr1=myrand(&(st->site_prn));
	      xr1 =  (log((double)(xr1+ 1.e-10)));
	      
	      xr2=myrand(&(st->site_prn));
	      xr2 =  (log((double)(xr2+ 1.e-10)));
	      
	      xr3=myrand(&(st->site_prn));
	      xr4=myrand(&(st->site_prn));

	      xr3=cos((double)pi2*xr3);

	      d = -(xr2 + xr1*xr3*xr3)/al;
	      if((1.00 - 0.5*d) > xr4*xr4) test = 1;
	    }
	    if(this_node == 0 && test !=1)
	      printf("site  took 20 kp hits\n");
	  } /* endif nacd */


/*------- creutz algorithm --------*/
	  if(nacd == 0 && al <= 2.0) {
printf("Mike's method\n");
	    cr++;
	    xl=exp((double)(-2.0*al));
	    xd= 1.0 - xl;
	    test=0;
	    for(k=0;k<20 && test == 0  ;k++) {
/* get two random numbers */
	      xr1=myrand(&(st->site_prn));
	      xr2=myrand(&(st->site_prn));
	      
	      r = xl + xd*xr1; 
	      a0 = 1.00 + log((double)r)/al;
	      if((1.0 -a0*a0) > xr2*xr2) test = 1;
	    }
	    d = 1.0 - a0;
	    if(this_node == 0 && test !=1) 
	      printf("site  took 20 creutz hits\n");
	  } /* endif nacd */


/*============ generate full su(2) matrix and update link matrix*/

/* find a0  = 1 - d*/
	  a0 = 1.0 - d;
/* compute r */
	  r2 = 1.0 - a0*a0;
	  r2 = fabs((double)r2);
	  r = sqrt((double)r2);

/* compute a3 */
	  a3=(2.0*myrand(&(st->site_prn)) - 1.0)*r;

/* compute a1 and a2 */
	  rho = r2 - a3*a3;
	  rho = fabs((double)rho);
	  rho= sqrt((double)rho);

/*xr2 is a random number between 0 and 2*pi */
	  xr2=pi2*myrand(&(st->site_prn));
	  
	  a1= rho*cos((double)xr2);
	  a2= rho*sin((double)xr2);

/*============ now do the updating.  h = a*action^dagger, new u = h*u */

	  a = make_su2_matrix(a0,a1,a2,a3);
printf("a = \n");
dump_su2_mat(&a);

	  mult_su2_na(&a,&action,&h);

/* update the link */
	  mult_su2_nn(&h,&(st->link[dir]),&a);
	  su2mat_copy(&a,&(st->link[dir]));
      
	  } /*   st */

/* diagnostics 
avekp=(float)kp / (float)(nx*ny*nz*nt/2);
avecr=(float)cr / (float)(nx*ny*nz*nt/2);
if(this_node ==0)
printf(" ave kp steps = %e, ave creutz steps = %e\n",
(double)avekp,(double)avecr);
*/
      } /*  direction */
    }} /* parity and NumTrj */
} /* monte */




