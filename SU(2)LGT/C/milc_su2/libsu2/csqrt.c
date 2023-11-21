/* Subroutines for operations on complex numbers */
/* complex square root */
#include <math.h>
#include "complex.h"

complex csqrt(z) complex *z; {
complex c;
float theta,r;
    r = sqrt(hypot(z->real,z->imag));
    theta = 0.5*atan2(z->imag,z->real);
    c = ce_itheta(theta);
    c.real *=r; c.imag *= r;
    return(c);
}
