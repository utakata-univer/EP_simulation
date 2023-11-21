/* Subroutines for operations on complex numbers */
/* Divide two complex numbers */
#include "complex.h"

complex cdiv(a,b) complex *a; complex *b;  {
    complex c;
    float scale;
    scale = 1.0/((*b).real*(*b).real+(*b).imag*(*b).imag);
    c.real = scale*((*a).real*(*b).real + (*a).imag*(*b).imag);
    c.imag = scale*((*a).imag*(*b).real - (*a).real*(*b).imag);
    return(c);
}
