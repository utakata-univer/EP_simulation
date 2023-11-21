/* Subroutines for operations on complex numbers */
/* double complex logarithm */
#include <math.h>
#include "complex.h"

double_complex dclog(a) double_complex *a;  {
    double_complex c;
    c.real = 0.5*(double)log((double)((*a).real*(*a).real+(*a).imag*(*a).imag));
    c.imag = (double)atan2( (double)(*a).imag, (double)(*a).real );
    return(c);
}
