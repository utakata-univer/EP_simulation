/*********************** reunitarize.c ***************************/
/* MIMD version 3 */

/* reunitarize SU(2) link matrices.
   Note: reunit_su2(mat) is in su2.a 
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "globaldefs.h"
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

void reunitarize() {
  register int i,dir;
  register site *s;
  register su2_matrix *mat;

  FORALLSITES(i,s)
    for(dir=XUP; dir<=TUP; dir++ ){
      mat = (su2_matrix *)&(s->link[dir]);
      reunit_su2( mat );
    }
}
