/* Subroutines for operations on complex numbers */
/* make a double complex number from two double precision reals */
#include "complex.h"

double_complex dcmplx(x,y) double x,y;  {
    double_complex c;
    c.real = x; c.imag = y;
    return(c);
}
