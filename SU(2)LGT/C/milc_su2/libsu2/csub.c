/* Subroutines for operations on complex numbers */
/* complex subtract */
#include "complex.h"

complex csub(a,b) complex *a; complex *b;  {
    complex c;
    c.real = (*a).real - (*b).real;
    c.imag = (*a).imag - (*b).imag;
    return(c);
}
