/* Subroutines for operations on complex numbers */
/* complex conjugate */
#include "complex.h"

complex conjg(a) complex *a;  {
    complex c;
    c.real = (*a).real;
    c.imag = -(*a).imag;
    return(c);
}
