#ifndef LATTICE_H
#define LATTICE_H

#define PI 3.14159265358979323846
#define L 10
#define D 4
#define V (L*L*L*L)

#ifdef CONTROL 
#define EXTERN 
#undef CONTROL
#else
#define EXTERN extern
#endif
 
EXTERN double a_0[V][4],
              a_1[V][4],
	      a_2[V][4],
	      a_3[V][4];
EXTERN int    hop[V][2*D];

#endif
