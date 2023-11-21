#ifndef SU2_H
#define SU2_H
#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"

void gauss_rand(int n, double* rand);
void hopping(int h[V][2*D]);

void print_meas();
void measure();
void init_measure();


/*data structure to store all the parameters of the algorithm*/
typedef struct {
               int ntherm ;     /*number of thermalization steps*/
	       int nsweep ;     /*number of sweep after thermalization*/
	       int naccu;       /*binsize for printing out the measurements*/
              } heatbath_params_t;

/*data structure to store all the parameters of the action*/
typedef struct {
              double beta;
              } act_params_t;


/*
double hmc(act_params_t *apara, hmc_params_t *hpara);
*/
double plaquette(void);
double heatbath(act_params_t *apara, heatbath_params_t *hpara);
#endif
