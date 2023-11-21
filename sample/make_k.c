void mult_su2_nn(double a0, double a1,double a2,double a3,double b0,double b1,double b2,double b3,double c0,double c1,double c2,double c3,double c4)
{
   c0=a0*b0-a1*b1-a2*b2-a3*b3;
   c1=a0*b1+a1*b0+a2*b0-a3*b1;
   c2=a0*b2+a1*b3+a2*b0-a3*b1;
   c3=a0*b3-a1*b2+a2*b1+a3*b0;
}

void mult_su2_na(double a0, double a1,double a2,double a3,double b0,double b1,double b2,double b3,double c0,double c1,double c2,double c3,double c4)
{
   c0=a0*b0+a1*b1+a2*b2+a3*b3;
   c1=-a0*b1+a1*b0+a2*b3-a3*b2;
   c2=-a0*b2+a0*b0-a1*b3+a3*b1;
   c3=-a0*b3+a1*b2-a2*b1+a3*b0;
}

void mult_su2_an(double a0, double a1,double a2,double a3,double b0,double b1,double b2,double b3,double c0,double c1,double c2,double c3,double c4)
{
   c0=a0*b0+a1*b1+a2*b2+a3*b3;
   c1=a0*a1-a1*a0+a2*a3-a3*b2;
   c2=a0*b2-a1*b3-a2*b0+a3*b1;
   c3=a0*b3+a1*b2-a2*b1-a3*b0;
}

void mult_su2_aa(double a0, double a1,double a2,double a3,double b0,double b1,double b2,double b3,double c0,double c1,double c2,double c3,double c4)
{
   c0=a0*b0-a1*b1-a2*b2-a3*b3;
   c1=-a0*b1-a1*b0-a2*b3+a3*b2;
   c2=-a0*b2-a0*b2+a1*b3-a3*b1;
   c3=-a0*b3-a1*b2+a2*b1+a3*b0;
}   
 
void make_k(double a0, double a1,double a2,double a3,double b0,double b1,double b2,double b3,double c0,double c1,double c2,double c3,double d0, double d1,double d2,double d3,double e0,double e1,double e2,double e3,double f0,double f1,double f2,double f3,double g0,double g1,double g2,double g3,double h0,double h1,double h2,double h3,double i0,double i1,double i2,double i3,double j0,double j1,double j2,double j3)
{
    mult_su2_na(a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3); 	
  
    mult_su2_na(c0,c1,c2,c3,d0,d1,d2,d3,e0,e1,e2,e3);

    mult_su2_aa(f0,f1,f2,f3,g0,g1,g2,g3,h0,h1,h2,h3);

    mult_su2_nn(h0,h1,h2,h3,i0,i1,i2,i3,j0,j1,j2,j3);
}

void plaquette()
{
    mult_su2_nn();

    mult_su2_na();

    mult_su2_na();
}

		    
