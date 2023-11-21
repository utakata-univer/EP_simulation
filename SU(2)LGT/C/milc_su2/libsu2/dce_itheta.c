/* Subroutines for operations on complex numbers */
/* double complex exp( i*theta ) */
#include <math.h>
#include "complex.h"

double_complex dce_itheta(theta) double theta;  {
    double_complex c;
    c.real = (double)cos( (double)theta );
    c.imag = (double)sin( (double)theta );
    /* there must be a more efficient way */
    return( c );
}
