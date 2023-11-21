#include "complex.h"
#include "globaldefs.h"
#include "su2.h"
#include <stdio.h>
#include <stdlib.h>

double sqrt(double);

int main() {
   su2_matrix a, b, c;
   su2_vector u, v, w;
   float scale;
   int i;

   srand(1234);
/*
   u.c[0].real = 1;
   u.c[0].imag = 1;
   u.c[1].real = 2;
   u.c[1].imag = 3;
*/
    u.c[0] = cmplx(1., 1.);
    u.c[1] = cmplx(1., 1.);

   
   dump_su2_vec(&u);

   a.e[0] = 2.;
   a.e[1] = 3.;
   a.e[2] = 7.;
   a.e[3] =13;

   scale = sqrt( a.e[0]*a.e[0]+a.e[1]*a.e[1]+a.e[2]*a.e[2]+a.e[3]*a.e[3]);
   for(i=0;i<=3;i++) a.e[i] /= scale;

   dump_su2_mat(&a);

   mult_su2_mat_vec(&a, &u, &w);
   printf("%e %e\n", magsq_su2vec(&u), magsq_su2vec(&w));  

/*
  for(i=1;i<=3;i++)b.e[i] = - a.e[i];
  
   mult_su2_nn(&a,&b,&c);
   dump_su2_mat(&c);
   mult_su2_an(&a,&a,&c);
   dump_su2_mat(&c);
   mult_su2_na(&a,&a,&c);
   dump_su2_mat(&c);
   mult_su2_aa(&a,&b,&c);
   dump_su2_mat(&c);
   mult_su2_choose(&a,1,&a,-1,&c);
   dump_su2_mat(&c);
*/
}
