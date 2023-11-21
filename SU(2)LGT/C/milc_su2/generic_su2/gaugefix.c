/************************** gaugefix.c *******************************/
/* Fix Coulomb or Lorentz gauge by doing successive SU(2) gauge hits */
/* MIMD version 3 */
/* C. DeTar 10-22-90 */

/* Prototype...
   void gaugefix(gauge_dir,relax_boost,max_gauge_iter,gauge_fix_tol)
   int gauge_dir,max_gauge_iter;
   float relax_boost,gauge_fix_tol;


   gauge_dir     specifies the direction of the "time"-like hyperplane
                 for the purposes of defining Coulomb or Lorentz gauge 
      TUP    for evaluating propagators in the time-like direction
      ZUP    for screening lengths.
      8      for Lorentz gauge
      relax_boost	Overrelaxation parameter 
      max_gauge_iter	Maximum number of iterations 
      gauge_fix_tol	Stop if change is less than this 
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include LATDEF
#include <comdefs.h>

/* Generic definitions - could be useful elsewhere */
#define FORALLUPDIRBUT(direction,dir) \
   for(dir=XUP; dir<= TUP; dir++)if(dir != direction)
#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. Nulls EVENANDODD*/
/*    CDIF(a,b)         a -= b						      */
								/*  a -= b    */
#define CDIF(a,b) { (a).real -= (b).real; (a).imag -= (b).imag; }

typedef struct { complex e[2][2]; } su2_matrix;

void mult_su2_mat_vec_elem_n(u,x0,x1)
su2_matrix *u;
complex *x0, *x1;
{
  /* Multiplies the complex column spinor (x0, x1) by the SU(2) matrix u */
  /* and puts the result in (x0,x1).  */
  /* Thus x <- u * x          */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL(u->e[0][0], t0, z0);
  CMUL(u->e[0][1], t1, z1);
  CADD(z0, z1, *x0);
  CMUL(u->e[1][0], t0, z0);
  CMUL(u->e[1][1], t1, z1);
  CADD(z0, z1, *x1);

} /* mult_su2_mat_vec_elem_n */

void mult_su2_mat_vec_elem_a(u,x0,x1)
su2_matrix *u;
complex *x0,*x1;
{
  /* Multiplies the complex row spinor (x0, x1) by the adjoint of the */
  /* SU(2) matrix u and puts the result in (x0,x1).  */
  /* Thus x <-  x * u-adj       */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL_J(t0, u->e[0][0], z0);
  CMUL_J(t1, u->e[0][1], z1);
  CADD(z0, z1, *x0);
  CMUL_J(t0, u->e[1][0], z0);
  CMUL_J(t1, u->e[1][1], z1);
  CADD(z0, z1, *x1);

} /* mult_su2_mat_vec_elem_a */

void dumpsu2(u) su2_matrix *u; {
  int i,j;
  for(i=0;i<2;i++){
    for(j=0;j<2;j++)printf("(%.2e,%.2e)\t",
			  (double)u->e[i][j].real,(double)u->e[i][j].imag);
    printf("\n");
  }
  printf("\n");
}

void left_su2_hit_n(u,p,q,link)
     su2_matrix *u;
     su3_matrix *link;
     int p,q;
{
  /* link <- u * link */
  /* The 0 row of the SU(2) matrix u matches row p of the SU(3) matrix */
  /* The 1 row of the SU(2) matrix u matches row q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  register int m;

  for (m = 0; m < 3; m++)
    mult_su2_mat_vec_elem_n(u, &(link->e[p][m]), &(link->e[q][m]));

} /* left_su2_hit_n */

void right_su2_hit_a(u,p,q,link)
     su2_matrix *u;
     su3_matrix *link;
     int p,q;
{
  /* link <-  link * u adj */
  /* The 0 column of u-adjoint matches column p of the SU(3) matrix */
  /* The 1 column of u-adjoint matches column q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  register int m;

  for (m = 0; m < 3; m++)
    mult_su2_mat_vec_elem_a(u, &(link->e[m][p]), &(link->e[m][q]));

} /*right_su2_hit_a */

void accum_gauge_hit(gauge_dir,parity)
     int gauge_dir;
     int parity;
{

/* Accumulates sums and differences of link matrices for determining optimum */
/* hit for gauge fixing */
/* Differences are kept in staple and the diagonal elements of the sums */
/* in tempvec  */

  register int j,k;
  register su3_matrix *m1,*m2;
  register int dir,i;
  register site *s;
  float *p2;

  /* Clear tempvec and staple */

  FORSOMEPARITY(i,s,parity)
    {
      for(j=0;j<3;j++)
	{
	  s->tempvec.c[j] = cmplx(0.,0.);
	  for(k=0;k<3;k++)
	      s->staple.e[j][k] = cmplx(0.,0.);
	}
    }
  
  /* Subtract upward link contributions */

  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  /* Upward link matrix */
	  m1 = &(s->link[dir]);
          sub_su3_matrix(&(s->staple), m1, &(s->staple));
          for(j=0;j<3;j++)CSUM(s->tempvec.c[j], m1->e[j][j]);
	}
    }

  /* Add downward link contributions */

  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  /* Downward link matrix */
	  m2 = (su3_matrix *)gen_pt[dir][i];
	  p2 = (float *)gen_pt[dir+4][i];

	  add_su3_matrix(&(s->staple), m2, &(s->staple));

	  /* Add diagonal elements to tempvec  */
	  for(j=0;j<3;j++)CSUM(s->tempvec.c[j], m2->e[j][j]);
	}
    }
} /* accum_gauge_hit */


void do_hit(gauge_dir, parity, p, q, relax_boost)
int gauge_dir, parity, p, q;
float relax_boost;
{
  /* Do optimum SU(2) gauge hit for p, q subspace */

  float a0,a1,a2,a3,asq,a0sq,x,r,xdr;
  register int dir,i;
  register site *s;
  register int m;
  su2_matrix u;
  register su3_matrix htemp;

  /* Accumulate sums for determining optimum gauge hit */

  accum_gauge_hit(gauge_dir,parity);

  FORSOMEPARITY(i,s,parity)
    {
      /* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
      /* The locally optimum unnormalized components a0, aj are determined */
      /* from the current link in direction dir and the link downlink */
      /* in the same direction on the neighbor in the direction opposite dir */
      /* The expression is */
      /* a0 = Sum dir Tr Re 1       * (downlink dir + link dir) */
      /* aj = Sum dir Tr Im sigma j * (downlink dir - link dir)  j = 1,2, 3 */
      /*   where 1, sigma j are unit and Pauli matrices on the p,q subspace */
      
      a0 =  s->tempvec.c[p].real + s->tempvec.c[q].real;
      a1 =  s->staple.e[q][p].imag + s->staple.e[p][q].imag;
      a2 = -s->staple.e[q][p].real + s->staple.e[p][q].real;
      a3 =  s->staple.e[p][p].imag - s->staple.e[q][q].imag;
      
      /* Over-relaxation boost */

      /* This algorithm is designed to give little change for large |a| */
      /* and to scale up the gauge transformation by a factor of relax_boost*/
      /* for small |a| */

      asq = a1*a1 + a2*a2 + a3*a3;
      a0sq = a0*a0;
      x = (relax_boost*a0sq + asq)/(a0sq + asq);
      r = sqrt((double)(a0sq + x*x*asq));
      xdr = x/r;
      /* Normalize and boost */
      a0 = a0/r; a1 = a1*xdr; a2 = a2*xdr; a3 = a3*xdr;

      /* Elements of SU(2) matrix */

      u.e[0][0] = cmplx( a0, a3);
      u.e[0][1] = cmplx( a2, a1);
      u.e[1][0] = cmplx(-a2, a1);
      u.e[1][1] = cmplx( a0,-a3);

      
      /* Do SU(2) hit on all upward links */

      FORALLUPDIR(dir)
	left_su2_hit_n(&u,p,q,&(s->link[dir]));
      
      /* Do SU(2) hit on all downward links */
      
      FORALLUPDIR(dir)
	right_su2_hit_a(&u,p,q,(su3_matrix *)gen_pt[dir][i]);
      
    }
  
  /* Exit with modified downward links left in communications buffer */
} /* do_hit */

float get_gauge_fix_action(gauge_dir,parity)
int gauge_dir;
int parity;
{
  /* Adds up the gauge fixing action for sites of given parity */
  /* Returns average over these sites */
  /* The average is normalized to a maximum of 1 when all */
  /* links are unit matrices */

  register int dir,i,ndir;
  register site *s;
  register su3_matrix *m1, *m2;
  float *p2;
  float gauge_fix_action;
  complex trace;

  gauge_fix_action = 0.0;
  
  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  m1 = &(s->link[dir]);
	  m2 = (su3_matrix *)gen_pt[dir][i];
	  p2 = (float *)gen_pt[dir+4][i];

	  trace = trace_su3(m1);
	  gauge_fix_action += trace.real;

	  trace = trace_su3(m2); 
 	  gauge_fix_action += trace.real;
	}
    }

  /* Count number of terms to average */
  ndir = 0; FORALLUPDIRBUT(gauge_dir,dir)ndir++;
  
  /* Sum over all sites of this parity */
  g_floatsum( &gauge_fix_action);
  
  /* Average is normalized to max of 1/2 on sites of one parity */
  return(gauge_fix_action /((float)(6*ndir*nx*ny*nz*nt)));
} /* get_gauge_fix_action */

void gaugefixstep(gauge_dir,av_gauge_fix_action,relax_boost)
int gauge_dir;
float *av_gauge_fix_action;
float relax_boost;
{
  /* Carry out one iteration in the gauge-fixing process */

  int parity;
  msg_tag *mtag[8];
  float gauge_fix_action;
  register int dir,i;
  register site *s;

  register int j,k;
  su3_matrix *m1, *m2;
  
  /* Alternate parity to prevent interactions during gauge transformation */

  *av_gauge_fix_action = 0.;
  g_sync();
  fflush(stdout);
  
  for(parity = ODD; parity <= EVEN; parity++)
    {
      /* Start gathers of downward links */
      
      FORALLUPDIR(dir)
	{
	  mtag[dir] = start_gather( F_OFFSET(link[dir]), sizeof(su3_matrix),
			   OPP_DIR(dir), parity, gen_pt[dir] );
	  /* Just in case someone changes the direction definitions */
          if(dir+4 >= 8)
 	      {printf("GFIX: Dimension overrun in mtag\n");terminate(1);}
	}
      
      /* Wait for gathers */
      
      FORALLUPDIR(dir)
         {
	  wait_gather(mtag[dir]);
	 }

      /* Total gauge fixing action for sites of this parity: Before */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      /* Do optimum gauge hit on various subspaces */

      do_hit(gauge_dir,parity,0,1, relax_boost);
      do_hit(gauge_dir,parity,1,2, relax_boost);
      do_hit(gauge_dir,parity,2,0, relax_boost);

      /* Total gauge fixing action for sites of this parity: After */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);
      
      *av_gauge_fix_action += gauge_fix_action;

      /* Scatter downward link matrices by gathering to sites of */
      /* opposite parity */

      FORALLUPDIR(dir)
	{
	  /* Synchronize before scattering to be sure the new modified link */
	  /* matrices are all ready to be scattered and staple is not */
	  /* overwritten before it is used */
	  g_sync();

	  /* First copy modified link for this dir */
	  /* from comm buffer or node to staple */

	  FORSOMEPARITY(i,s,parity)
	    {
	      su3mat_copy((su3_matrix *)(gen_pt[dir][i]),&(s->staple));
	    }
	  
	  /* Now we are finished with gen_pt[dir] and gen_pt[dir+4] */
	  cleanup_gather(mtag[dir]);
      
	  /* Synchronize to make sure the previous copy happens before the */
	  /* subsequent gather below  */
	  g_sync();

	  /* Gather staple onto sites of opposite parity */
	  mtag[dir] = start_gather( F_OFFSET(staple), sizeof(su3_matrix),
			    dir, OPP_PAR(parity), gen_pt[dir] );
	  wait_gather(mtag[dir]);

         /* Copy modified matrices into proper location */

         FORSOMEPARITY(i,s,OPP_PAR(parity))
	      su3mat_copy((su3_matrix *)(gen_pt[dir][i]),&(s->link[dir]));

         cleanup_gather(mtag[dir]);
	}

    }
} /* gaugefix */

void gaugefix(gauge_dir,relax_boost,max_gauge_iter,gauge_fix_tol)
int gauge_dir,max_gauge_iter;
float relax_boost,gauge_fix_tol;
{
  int gauge_iter;
  float current_av, old_av, del_av;

  /* Do at most max_gauge_iter iterations, but stop after the second step if */
  /* the change in the avg gauge fixing action is smaller than gauge_fix_tol */

  for (gauge_iter=0; gauge_iter < max_gauge_iter; gauge_iter++)
    {
      gaugefixstep(gauge_dir,&current_av,relax_boost);

      if(gauge_iter != 0)
	{
	  del_av = current_av - old_av;
	  if (fabs(del_av) < gauge_fix_tol) break;
	}
      old_av = current_av;
    }
  if(this_node==0)
    printf("GFIX: Ended at step %d. Av gf action %.3e, delta %.3e\n",
	   gauge_iter,(double)current_av,(double)del_av);
}
