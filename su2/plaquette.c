#include "su2.h"

/*
static int d[4][4] = {

	{1,0,0,0},
	{0,1,0,0},
	{0,0,1,0},
	{0,0,0,1},
};
  

void kei(void)
{
   int nx,ny,nz,nt,mu,nu;	   
   
       for (nx=1;nx<L+1;nx++)
    {
       for (ny=1;ny<L+1;ny++)
    {
       for (nz=1;nz<L+1;nz++)
    {
       for (nt=1;nt<L+1;nt++)
    {
       for (mu=0;mu<4;mu++)
    {
       a_0[0][ny][nz][nt][mu]=a_0[L][ny][nz][nt][mu];
       a_1[0][ny][nz][nt][mu]=a_1[L][ny][nz][nt][mu];
       a_2[0][ny][nz][nt][mu]=a_2[L][ny][nz][nt][mu];
       a_3[0][ny][nz][nt][mu]=a_3[L][ny][nz][nt][mu];
       a_0[nx][0][nz][nt][mu]=a_0[nx][L][nz][nt][mu];
       a_1[nx][0][nz][nt][mu]=a_1[nx][L][nz][nt][mu];
       a_2[nx][0][nz][nt][mu]=a_2[nx][L][nz][nt][mu];
       a_3[nx][0][nz][nt][mu]=a_3[nx][L][nz][nt][mu];
       a_0[nx][ny][0][nt][mu]=a_0[nx][ny][L][nt][mu];
       a_1[nx][ny][0][nt][mu]=a_1[nx][ny][L][nt][mu];
       a_2[nx][ny][0][nt][mu]=a_2[nx][ny][L][nt][mu];
       a_3[nx][ny][0][nt][mu]=a_3[nx][ny][L][nt][mu];
       a_0[nx][ny][nz][0][mu]=a_0[nx][ny][nz][L][mu];
       a_1[nx][ny][nz][0][mu]=a_1[nx][ny][nz][L][mu];
       a_2[nx][ny][nz][0][mu]=a_2[nx][ny][nz][L][mu];
       a_3[nx][ny][nz][0][mu]=a_3[nx][ny][nz][L][mu];
       a_0[L+1][ny][nz][nt][mu]=a_0[1][ny][nz][nt][mu];
       a_1[L+1][ny][nz][nt][mu]=a_1[1][ny][nz][nt][mu];
       a_2[L+1][ny][nz][nt][mu]=a_2[1][ny][nz][nt][mu];
       a_3[L+1][ny][nz][nt][mu]=a_3[1][ny][nz][nt][mu];
       a_0[nx][L+1][nz][nt][mu]=a_0[nx][1][nz][nt][mu];
       a_1[nx][L+1][nz][nt][mu]=a_1[nx][1][nz][nt][mu];
       a_2[nx][L+1][nz][nt][mu]=a_2[nx][1][nz][nt][mu];
       a_3[nx][L+1][nz][nt][mu]=a_3[nx][1][nz][nt][mu];
       a_0[nx][ny][L+1][nt][mu]=a_0[nx][ny][1][nt][mu];
       a_1[nx][ny][L+1][nt][mu]=a_1[nx][ny][1][nt][mu];
       a_2[nx][ny][L+1][nt][mu]=a_2[nx][ny][1][nt][mu];
       a_3[nx][ny][L+1][nt][mu]=a_3[nx][ny][1][nt][mu];
       a_0[nx][ny][nz][L+1][mu]=a_0[nx][ny][nz][1][mu];
       a_1[nx][ny][nz][L+1][mu]=a_1[nx][ny][nz][1][mu];
       a_2[nx][ny][nz][L+1][mu]=a_2[nx][ny][nz][1][mu];
       a_3[nx][ny][nz][L+1][mu]=a_3[nx][ny][nz][1][mu];
     }

     }

     }

     }

     }
  

       for (nx=1;nx<L+1;nx++)
    {
       for (ny=1;ny<L+1;ny++)
    {
       for (nz=1;nz<L+1;nz++)
    {
       for (nt=1;nt<L+1;nt++)
    {
       for (mu=0;mu<4;mu++)    
    {   
       for (nu=0;nu<4;nu++)
    {	       

    p_0[nx][ny][nz][nt][mu][nu]=a_0[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_1[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_2[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_3[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];
   
    p_1[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_1[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_2[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-a_3[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];
   
    p_2[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-a_1[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_2[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_3[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];
   
    p_3[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_1[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-a_2[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_3[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];
   
    q_0[nx][ny][nz][nt][mu][nu]=p_0[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu];
 
    q_1[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]-p_3[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu];
    q_2[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]-p_1[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu];

    q_3[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]-p_2[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu];

    r_0[nx][ny][nz][nt][mu][nu]=a_0[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_1[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_2[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_3[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu];
       
    r_1[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]+a_1[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_2[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[3][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]+a_3[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu];
       
    r_2[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]+a_1[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_2[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_3[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu];
 
    r_3[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_1[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_2[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu]-a_3[nx+d[0][mu]-d[0][nu]][ny+d[1][mu]-d[1][nu]][nz+d[2][mu]-d[2][nu]][nt+d[3][mu]-d[3][nu]][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][mu];

    s_0[nx][ny][nz][nt][mu][nu]=r_0[nx][ny][nz][nt][mu][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_1[nx][ny][nz][nt][mu][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_2[nx][ny][nz][nt][mu][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_3[nx][ny][nz][nt][mu][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu];

    s_1[nx][ny][nz][nt][mu][nu]=r_0[nx][ny][nz][nt][mu][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_1[nx][ny][nz][nt][mu][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_2[nx][ny][nz][nt][mu][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_3[nx][ny][nz][nt][mu][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu];
   
    s_2[nx][ny][nz][nt][mu][nu]=r_0[nx][ny][nz][nt][mu][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_1[nx][ny][nz][nt][mu][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_2[nx][ny][nz][nt][mu][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_3[nx][ny][nz][nt][mu][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu];
 
    s_3[nx][ny][nz][nt][mu][nu]=r_0[nx][ny][nz][nt][mu][nu]*a_3[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]-r_1[nx][ny][nz][nt][mu][nu]*a_2[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_2[nx][ny][nz][nt][mu][nu]*a_1[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu]+r_3[nx][ny][nz][nt][mu][nu]*a_0[nx-d[0][nu]][ny-d[1][nu]][nz-d[2][nu]][nt-d[3][nu]][nu];

    }

    }

    }

    }

    }

    }

       for (nx=1;nx<L+1;nx++)
    {
       for (ny=1;ny<L+1;ny++)
    {
       for (nz=1;nz<L+1;nz++)
    {
       for (nt=1;nt<L+1;nt++)
    {
       for (mu=0;mu<4;mu++)
    {

	    k_0[nx][ny][nz][nt][mu]=q_0[nx][ny][nz][nt][mu][0]+s_0[nx][ny][nz][nt][mu][0]+q_0[nx][ny][nz][nt][mu][1]+s_0[nx][ny][nz][nt][mu][1]+q_0[nx][ny][nz][nt][mu][2]+s_0[nx][ny][nz][nt][mu][2]+q_0[nx][ny][nz][nt][mu][3]+s_0[nx][ny][nz][nt][mu][3]-q_0[nx][ny][nz][nt][mu][mu]-s_0[nx][ny][nz][nt][mu][mu];
 
	    k_1[nx][ny][nz][nt][mu]=q_1[nx][ny][nz][nt][mu][0]+s_1[nx][ny][nz][nt][mu][0]+q_1[nx][ny][nz][nt][mu][1]+s_1[nx][ny][nz][nt][mu][1]+q_1[nx][ny][nz][nt][mu][2]+s_1[nx][ny][nz][nt][mu][2]+q_1[nx][ny][nz][nt][mu][3]+s_1[nx][ny][nz][nt][mu][3]-q_1[nx][ny][nz][nt][mu][mu]-s_1[nx][ny][nz][nt][mu][mu];
 
	    k_2[nx][ny][nz][nt][mu]=q_2[nx][ny][nz][nt][mu][0]+s_2[nx][ny][nz][nt][mu][0]+q_2[nx][ny][nz][nt][mu][1]+s_2[nx][ny][nz][nt][mu][1]+q_2[nx][ny][nz][nt][mu][2]+s_2[nx][ny][nz][nt][mu][2]+q_2[nx][ny][nz][nt][mu][3]+s_2[nx][ny][nz][nt][mu][3]-q_2[nx][ny][nz][nt][mu][mu]-s_2[nx][ny][nz][nt][mu][mu];

	    k_3[nx][ny][nz][nt][mu]=q_3[nx][ny][nz][nt][mu][0]+s_3[nx][ny][nz][		    
*/
double plaquette(void)
{
   int i,mu,nu;	
   double P,p_0,p_1,p_2,p_3,q_0,q_1,q_2,q_3,r_0,r_1,r_2,r_3,s_0;

   P=0.0;   
       for (i=0;i<V;i++)
    {
       for (mu=0;mu<3;mu++)
    {
       for (nu=mu+1;nu<4;nu++)
    {
      p_0=a_0[i][mu];
      p_1=a_1[i][mu];
      p_2=a_2[i][mu];
      p_3=a_3[i][mu];
 
      q_0=p_0*a_0[hop[i][mu]][nu]-p_1*a_1[hop[i][mu]][nu]-p_2*a_2[hop[i][mu]][nu]-p_3*a_3[hop[i][mu]][nu];

      q_1=p_0*a_1[hop[i][mu]][nu]+p_1*a_0[hop[i][mu]][nu]-p_2*a_3[hop[i][mu]][nu]+p_3*a_2[hop[i][mu]][nu];

      q_2=p_0*a_2[hop[i][mu]][nu]+p_1*a_3[hop[i][mu]][nu]+p_2*a_0[hop[i][mu]][nu]-p_3*a_1[hop[i][mu]][nu];

      q_3=p_0*a_3[hop[i][mu]][nu]-p_1*a_2[hop[i][mu]][nu]+p_2*a_1[hop[i][mu]][nu]+p_3*a_0[hop[i][mu]][nu];
 
      r_0=q_0*a_0[hop[i][nu]][mu]+q_1*a_1[hop[i][nu]][mu]+q_2*a_2[hop[i][nu]][mu]+q_3*a_3[hop[i][nu]][mu];

      r_1=-q_0*a_1[hop[i][nu]][mu]+q_1*a_0[hop[i][nu]][mu]+q_2*a_3[hop[i][nu]][mu]-q_3*a_2[hop[i][nu]][mu];

      r_2=-q_0*a_2[hop[i][nu]][mu]-q_1*a_3[hop[i][nu]][mu]+q_2*a_0[hop[i][nu]][mu]+q_3*a_1[hop[i][nu]][mu];

      r_3=-q_0*a_3[hop[i][nu]][mu]+q_1*a_2[hop[i][nu]][mu]-q_2*a_1[hop[i][nu]][mu]+q_3*a_0[hop[i][nu]][mu];

      s_0=r_0*a_0[i][nu]+r_1*a_1[i][nu]+r_2*a_2[i][nu]+r_3*a_3[i][nu];

      P+=(1.0-s_0)/(6.0*V);

    }

    }

    }

    return P;

}    

/*

      p_0=p_0*a_0[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu]-p_1*a_1[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu]-p_2*a_2[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu]-p_3*a_3[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu];

      p_1=-p_0*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+p_1*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+p_2*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-p_3*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];
 
      p_2=-p_0*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-p_1*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+p_2*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+p_3*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];

      p_3[nx][ny][nz][nt][mu][nu]=-a_0[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_3[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_1[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_2[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]-a_2[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_1[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu]+a_3[nx+d[0][mu]][ny+d[1][mu]][nz+d[2][mu]][nt+d[3][mu]][nu]*a_0[nx+d[0][nu]][ny+d[1][nu]][nz+d[2][nu]][nt+d[3][nu]][mu];

     q_0[nx][ny][nz][nt][mu][nu]=p_0[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu];

    q_1[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]-p_3[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu];
    q_2[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]-p_1[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]+p_2[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu];

    q_3[nx][ny][nz][nt][mu][nu]=-p_0[nx][ny][nz][nt][mu][nu]*a_3[nx][ny][nz][nt][nu]+p_1[nx][ny][nz][nt][mu][nu]*a_2[nx][ny][nz][nt][nu]-p_2[nx][ny][nz][nt][mu][nu]*a_1[nx][ny][nz][nt][nu]+p_3[nx][ny][nz][nt][mu][nu]*a_0[nx][ny][nz][nt][nu];
      pla_0[nx][ny][nz][nt][mu][nu]=2*(a_0[nx][ny][nz][nt][mu]*q_0[nx][ny][nz][nt][mu][nu]-a_1[nx][ny][nz][nt][mu]*q_1[nx][ny][nz][nt][mu][nu]-a_2[nx][ny][nz][nt][mu]*q_2[nx][ny][nz][nt][mu][nu]-a_3[nx][ny][nz][nt][mu]*q_3[nx][ny][nz][nt][mu][nu]);
    }

    }

    }   
    
    }

    }

    }
    
    plaquette=0;

      for  (nx=1;nx<L+1;nx++)
    {
       for (ny=1;ny<L+1;ny++)
    {
       for (nz=1;nz<L+1;nz++)
    {
       for (nt=1;nt<L+1;nt++)
    {

      plaquette+=(pla_0[nx][ny][nz][nt][0][1]+pla_0[nx][ny][nz][nt][0][2]+pla_0[nx][ny][nz][nt][0][3]+pla_0[nx][ny][nz][nt][1][0]+pla_0[nx][ny][nz][nt][1][2]+pla_0[nx][ny][nz][nt][1][3]+pla_0[nx][ny][nz][nt][2][0]+pla_0[nx][ny][nz][nt][2][1]+pla_0[nx][ny][nz][nt][2][3]+pla_0[nx][ny][nz][nt][3][0]+pla_0[nx][ny][nz][nt][3][1]+pla_0[nx][ny][nz][nt][3][2])/(12*V);

    }

    }

    }

    }   

     
    return plaquette;
   
}

*/

