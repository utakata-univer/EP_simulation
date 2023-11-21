void triala_0(double a,double k,)
{
  ranlxd(&c,1);

  b=c*(1-exp(-2*beta*k))+exp(-2*beta*k);

  a=1+log(b)/(beta*k);

  delta=sqrt(1-a*a);
}


