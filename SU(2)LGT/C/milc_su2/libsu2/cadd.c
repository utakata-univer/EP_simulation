/* Subroutines for operations on complex numbers */
/* add two complex numbers */
#include "complex.h"

complex cadd(a,b) complex *a; complex *b;  {
    complex c;
    c.real = (*a).real + (*b).real;
    c.imag = (*a).imag + (*b).imag;
    return(c);
}
