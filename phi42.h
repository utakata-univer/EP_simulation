#ifndef PHI42_H
#define PHI42_H
#include "stdio.h"
#include "stdlib.h"
#include "lattice2.h"
#include "math.h"

void hopping(int h[V][2*D]);

typedef struct {
	       double tlength;
	       int nstep;
	       int ntherm ;
               int ntraj ;
              } hmc_params_t;


#endif
