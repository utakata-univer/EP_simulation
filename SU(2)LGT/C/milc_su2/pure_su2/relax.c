/************************** relax.c *******************************/
/* Microcanonical overrelaxation for SU(2) */
/* MIMD version 4 */
/* C. DeTar 18 Oct 1990 */
/* T. DeGrand March 1991 */
/* J. Hetrick and D. Toussaint April 1995 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>


/* Generic definitions */
#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)

void relax(int NumStp) {
  int NumTrj, parity;
  register int dir,i;
  register site *st;
  su2_matrix action, u, tmp;
  void dsdu_qhb();

  for(NumTrj=0; NumTrj<NumStp; NumTrj++)
    for(parity=ODD;parity<=EVEN;parity++)
      FORALLUPDIR(dir) {
/* compute the gauge force */
	dsdu_qhb(dir,parity);

/* Over-relax */
	FORSOMEPARITY(i,st,parity) {
	  reunit_su2(&(st->staple));
	  mult_su2_na(&(st->link[dir]), &(st->staple), &action);
	  mult_su2_an(&action, &(st->link[dir]), &tmp);
	  mult_su2_an(&action, &tmp, &(st->link[dir]));
	}

      }
   
} /* relax */

