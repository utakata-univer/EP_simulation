/* Subroutines for operations on complex numbers */
/* multiply two complex numbers */
#include "complex.h"

complex cmul(a,b) complex *a; complex *b;  {
    complex c;
    c.real = (*a).real * (*b).real - (*a).imag * (*b).imag;
    c.imag = (*a).imag * (*b).real + (*a).real * (*b).imag;
    return(c);
}
