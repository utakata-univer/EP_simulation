/******************  rand_ahmat.c  (in su2.a) ***************************
*									*
* void random_anti_hermitian( su2_anti_hermitmat *mat_antihermit, 
                                                       passthru *prn_pt)*
* Creates gaussian random anti-hermitian matrices			*
* Normalization is < e[a]*e[a] > = 1/2                              	*
* The argument "prn_pt" is a pointer to be passed to gaussian_rand_no() *
*/

typedef char passthru;
#include <math.h>
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

void random_anti_hermitian(su2_anti_hermitmat *mat_antihermit, 
			   void *prn_pt) {

	mat_antihermit->e[0] = gaussian_rand_no(prn_pt);
	mat_antihermit->e[1] = gaussian_rand_no(prn_pt);
	mat_antihermit->e[2] = gaussian_rand_no(prn_pt);

}/*random_anti_hermitian_*/
