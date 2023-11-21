/* Subroutines for operations on complex numbers */
/* double complex exponential */
#include <math.h>
#include "complex.h"

double_complex dcexp(a) double_complex *a;  {
    double_complex c;
    double mag;
    mag = (double)exp( (double)(*a).real );
    c.real = mag*(double)cos( (double)(*a).imag );
    c.imag = mag*(double)sin( (double)(*a).imag );
    return(c);
}
