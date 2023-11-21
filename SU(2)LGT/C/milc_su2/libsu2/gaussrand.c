/*****************  gaussrand.c  (in su2.a) *****************************
*									*
*  float gaussian_ran_no( passthru *prn_pt )				*
*  Gaussian distributed random number					*
*  Probability distribution exp( -x*x ), so < x^2 > = 1/2		*
*  This requires a random number generator named "myrand()", returning	*
*  a float uniformly distributed between zero and one. The argument of	*
*  this routine is a pointer to be passed to myrand(). 			*
*/
typedef char passthru;
#include <math.h>
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

float gaussian_rand_no(void *prn_pt) {
   float myrand();
   static int iset=0;
   static float gset;
   float fac,r,v1,v2;

    if  (iset == 0) {
	do {
	    v1=2.0*myrand(prn_pt)-1.0;
	    v2=2.0*myrand(prn_pt)-1.0;
	    r=v1*v1+v2*v2;
	} while (r >= 1.0);
	fac=sqrt( -log((double)r)/(double)r);
	gset=v1*fac;
	iset=1;
	return v2*fac;
    } else {
	iset=0;
	return gset;
    }
}

