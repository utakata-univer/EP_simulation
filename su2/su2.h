#ifndef SU2_H
#define SU2_H
#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"

void initialize();
void hopping(int h[V][2*D]);
double plaquette();

typedef struct {
	       double beta;
	       int ntraj;	       
	      } hb_params_t;

double heatbath(hb_params_t *hpara);

#endif
