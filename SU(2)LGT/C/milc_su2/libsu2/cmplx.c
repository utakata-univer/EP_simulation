/* Subroutines for operations on complex numbers */
/* make a complex number from two real numbers */
#include "complex.h"

complex cmplx(x,y) float x,y;  {
    complex c;
    c.real = x; c.imag = y;
    return(c);
}
