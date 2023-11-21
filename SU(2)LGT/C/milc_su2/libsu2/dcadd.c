/* Subroutines for operations on complex numbers */
/* double complex add */
#include "complex.h"

double_complex dcadd(a,b) double_complex *a; double_complex *b;  {
    double_complex c;
    c.real = (*a).real + (*b).real;
    c.imag = (*a).imag + (*b).imag;
    return(c);
}
