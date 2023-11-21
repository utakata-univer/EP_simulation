/* Subroutines for operations on complex numbers */
/* double complex multiply */
#include "complex.h"

double_complex dcmul(a,b) double_complex *a; double_complex *b;  {
    double_complex c;
    c.real = (*a).real * (*b).real - (*a).imag * (*b).imag;
    c.imag = (*a).imag * (*b).real + (*a).real * (*b).imag;
    return(c);
}
