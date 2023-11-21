/* Subroutines for operations on complex numbers */
/* double precision complex conjugate */
#include "complex.h"

double_complex dconjg(a) double_complex *a;  {
    double_complex c;
    c.real = (*a).real;
    c.imag = -(*a).imag;
    return(c);
}
