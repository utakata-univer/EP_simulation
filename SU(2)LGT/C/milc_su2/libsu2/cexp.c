/* Subroutines for operations on complex numbers */
/* complex exponential */
#include <math.h>
#include "complex.h"

complex cexp(a) complex *a;  {
    complex c;
    float mag;
    mag = (float)exp( (double)(*a).real );
    c.real = mag*(float)cos( (double)(*a).imag );
    c.imag = mag*(float)sin( (double)(*a).imag );
    return(c);
}
