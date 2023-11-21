/* Subroutines for operations on complex numbers */
/* exp( i*theta ) */
#include <math.h>
#include "complex.h"

complex ce_itheta(theta) float theta;  {
    complex c;
    c.real = (float)cos( (double)theta );
    c.imag = (float)sin( (double)theta );
    /* there must be a more efficient way */
    return( c );
}
