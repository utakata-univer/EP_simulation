#ifndef LATTICE_H
#define LATTICE_H
/* Dimension of the lattice */
#define D 4
/* spatial extend of the lattice */
#define L 8
/* lattice volume, needs to be adjusted according to number of dimensions*/
#define V (L*L*L*L)

#ifdef CONTROL 
#define EXTERN 
#undef CONTROL
#else
#define EXTERN extern
#endif

EXTERN double a[V][D];
EXTERN int    hop[V][2*D];
#endif
