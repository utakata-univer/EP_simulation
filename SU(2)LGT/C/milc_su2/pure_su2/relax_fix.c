/************************** relax.c *******************************/
/* Microcanonical overrelaxation by doing successive SU(2) gauge hits */
/* MIMD version 3 */
/* T. DeGrand March 1991 */


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include "lattice.h"
#include <comdefs.h>

#define Nc 3

/* Generic definitions - could be useful elsewhere */
#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. Nulls EVENANDODD*/

typedef struct { complex e[2][2]; } su2_matrix;

/* Codes for interpreting selection of gauge fixing options */


/*  define out routines mult_su2_mat_vec_elem_n(u,x0,x1),  */
/*   dumpsu2(u), and  left_su2_hit_n(u,p,q,link)            */
/*   --- these routines also appear in gaugefix.c           */
/*   --CWB 4/1/91                                           */
/*                                                          */
 
#ifdef DONTCOMMENTOUT


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

#endif

/*  define out routines mult_su2_mat_vec_elem_n(u,x0,x1),  */
/*   dumpsu2(u), and  left_su2_hit_n(u,p,q,link)            */
/*   --- these routines also appear in gaugefix.c           */
/*   --CWB 4/1/91                                           */
/*                                                          */




void relax(NumStp)
int NumStp;
{
  /* Do overrelaxation by SU(2) subgroups */
int NumTrj,Nhit, index1, ina, inb,ii;
int parity;
float a0,a1,a2,a3,asq,r;
register int dir,i;
register site *st;
int first,count;
su3_matrix action;  su2_matrix u;
void dsdu_qhb();

Nhit = 3;
	for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++)
	for(parity=ODD;parity<=EVEN;parity++)
	{
	FORALLUPDIR(dir)
		{
              /* compute the gauge force */
		dsdu_qhb(dir,parity);
              /* now for the overrelaxed updating */
		for(index1=0;index1<Nhit;index1++)
		{
                      /*  pick out an SU(2) subgroup */
			ina=(index1+1) % Nc;
			inb=(index1+2) % Nc;
			if(ina > inb){ ii=ina; ina=inb; inb=ii;}


		FORSOMEPARITY(i,st,parity){
			mult_su3_na( &(st->link[dir]), &(st->staple), &action );

/*decompose the action into SU(2) subgroups using Pauli matrix expansion */
/* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
			a0 =  action.e[ina][ina].real + action.e[inb][inb].real;
			a3 =  action.e[ina][ina].imag - action.e[inb][inb].imag;
			a1 =  action.e[ina][inb].imag + action.e[inb][ina].imag;
			a2 =  action.e[ina][inb].real - action.e[inb][ina].real;
      


      /* Normalize and complex conjugate u */

			asq = a0*a0 + a1*a1 + a2*a2 + a3*a3;
			r = sqrt((double)asq );
			a0 = a0/r; a1 = -a1/r; a2 = -a2/r; a3 = -a3/r;
      /* Elements of SU(2) matrix */

			u.e[0][0] = cmplx( a0, a3);
			u.e[0][1] = cmplx( a2, a1);
			u.e[1][0] = cmplx(-a2, a1);
			u.e[1][1] = cmplx( a0,-a3);
    
      /* Do SU(2) hit on all links twice (to overrelax)  */

			left_su2_hit_n(&u,ina,inb,&(st->link[dir]));
			left_su2_hit_n(&u,ina,inb,&(st->link[dir])); 

			} /*   st */
		} /*  hits */
		} /*  direction */
	} /* parity, NumTrj*/
 
} /* relax */

