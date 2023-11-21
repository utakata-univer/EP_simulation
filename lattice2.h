#ifndef LATTICE2_H
#define LATTICE2_H

#define D 3
#define L 8
#define V (L*L*L)

#ifdef CONTROL
#define EXTERN
#undef CONTROL
#else
#define EXTERN extern
#endif

EXTERN double phi[V];
EXTERN int    hop[V][2*D];

#endif

