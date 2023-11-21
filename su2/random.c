#include "su2.h"

static double beta;
/*
static int ntraj;

static int d[4][4] = {

        {1,0,0,0},
        {0,1,0,0},
        {0,0,1,0},
        {0,0,0,1},
};
*/
void initialize(void)
{
    int i,mu;
    
       for (i=0;i<V;i++){
       for (mu=0;mu<4;mu++){    
	    a_0[i][mu]=1.0;
            a_1[i][mu]=0.0;
            a_2[i][mu]=0.0;
            a_3[i][mu]=0.0;
                           }
                           }
}


static void update_heatbath(void)
{
    int i,mu,nu;
    double b,c,r,x,y,theta,phi,delta,k;                                         
    double r_0,r_1,r_2,r_3,k_0,k_1,k_2,k_3,l_0,l_1,l_2,l_3,m_0,m_1,m_2,m_3,n_0,n_1,n_2,n_3,o_0,o_1,o_2,o_3,p_0,p_1,p_2,p_3,q_0,q_1,q_2,q_3;	    

  
     for (i=0;i<V;i++)
    {
     for (mu=0;mu<4;mu++)
    {     
  
    q_0=0.0;
    q_1=0.0;
    q_2=0.0;
    q_3=0.0;

    for (nu=0;nu<4;nu++)
    {	   if(nu != mu){
	    
     k_0=a_0[hop[i][mu]][nu];
     k_1=a_1[hop[i][mu]][nu];
     k_2=a_2[hop[i][mu]][nu];
     k_3=a_3[hop[i][mu]][nu];
 
     l_0=k_0*a_0[hop[i][nu]][mu]+k_1*a_1[hop[i][nu]][mu]+k_2*a_2[hop[i][nu]][mu]+k_3*a_3[hop[i][nu]][mu];
     l_1=-k_0*a_1[hop[i][nu]][mu]+k_1*a_0[hop[i][nu]][mu]+k_2*a_3[hop[i][nu]][mu]-k_3*a_2[hop[i][nu]][mu];
     l_2=-k_0*a_2[hop[i][nu]][mu]-k_1*a_3[hop[i][nu]][mu]+k_2*a_0[hop[i][nu]][mu]+k_3*a_1[hop[i][nu]][mu];
     l_3=-k_0*a_3[hop[i][nu]][mu]+k_1*a_2[hop[i][nu]][mu]-k_2*a_1[hop[i][nu]][mu]+k_3*a_0[hop[i][nu]][mu]; 

     m_0=l_0*a_0[i][nu]+l_1*a_1[i][nu]+l_2*a_2[i][nu]+l_3*a_3[i][nu]; 
     m_1=-l_0*a_1[i][nu]+l_1*a_0[i][nu]+l_2*a_3[i][nu]-l_3*a_2[i][nu];
     m_2=-l_0*a_2[i][nu]-l_1*a_3[i][nu]+l_2*a_0[i][nu]+l_3*a_1[i][nu];
     m_3=-l_0*a_3[i][nu]+l_1*a_2[i][nu]-l_2*a_1[i][nu]+l_3*a_0[i][nu];

     q_0+=m_0;
     q_1+=m_1;
     q_2+=m_2;
     q_3+=m_3;
     
     n_0=a_0[hop[hop[i][mu]][D+nu]][nu];
     n_1=-a_1[hop[hop[i][mu]][D+nu]][nu];
     n_2=-a_2[hop[hop[i][mu]][D+nu]][nu];
     n_3=-a_3[hop[hop[i][mu]][D+nu]][nu];

     o_0=n_0*a_0[hop[i][D+nu]][mu]+n_1*a_1[hop[i][D+nu]][mu]+n_2*a_2[hop[i][D+nu]][mu]+n_3*a_3[hop[i][D+nu]][mu];
     o_1=-n_0*a_1[hop[i][D+nu]][mu]+n_1*a_0[hop[i][D+nu]][mu]+n_2*a_3[hop[i][D+nu]][mu]-n_3*a_2[hop[i][D+nu]][mu];
     o_2=-n_0*a_2[hop[i][D+nu]][mu]-n_1*a_3[hop[i][D+nu]][mu]+n_2*a_0[hop[i][D+nu]][mu]+n_3*a_1[hop[i][D+nu]][mu];
     o_3=-n_0*a_3[hop[i][D+nu]][mu]+n_1*a_2[hop[i][D+nu]][mu]-n_2*a_1[hop[i][D+nu]][mu]+n_3*a_0[hop[i][D+nu]][mu];

     p_0=o_0*a_0[hop[i][D+nu]][nu]-o_1*a_1[hop[i][D+nu]][nu]-o_2*a_2[hop[i][D+nu]][nu]-o_3*a_3[hop[i][D+nu]][nu];
     p_1=o_0*a_1[hop[i][D+nu]][nu]+o_1*a_0[hop[i][D+nu]][nu]-o_2*a_3[hop[i][D+nu]][nu]+o_3*a_2[hop[i][D+nu]][nu];
     p_2=o_0*a_2[hop[i][D+nu]][nu]+o_1*a_3[hop[i][D+nu]][nu]+o_2*a_0[hop[i][D+nu]][nu]-o_3*a_1[hop[i][D+nu]][nu];
     p_3=o_0*a_3[hop[i][D+nu]][nu]-o_1*a_2[hop[i][D+nu]][nu]+o_2*a_1[hop[i][D+nu]][nu]+o_3*a_0[hop[i][D+nu]][nu];

     q_0+=p_0;
     q_1+=p_1;
     q_2+=p_2;
     q_3+=p_3;
    }
    }
    k=sqrt(q_0*q_0+q_1*q_1+q_2*q_2+q_3*q_3);

    do {

    ranlxd(&b,1);

    c=(1.0-exp(-2.0*beta*k))*b+exp(-2.0*beta*k);

    r_0=1.0+log(c)/(beta*k);

    delta=sqrt(1.0-r_0*r_0);

    ranlxd(&r,1);

    } while (1.0-delta<r);
    
    ranlxd(&x,1);
    ranlxd(&y,1);

    theta=2.0*x-1.0;
    phi=2.0*PI*y;

    r_1=sqrt(1.0-a_0[i][mu]*a_0[i][mu])*sqrt(1.0-theta*theta)*cos(phi);
    r_2=sqrt(1.0-a_0[i][mu]*a_0[i][mu])*sqrt(1.0-theta*theta)*sin(phi);
    r_3=sqrt(1.0-a_0[i][mu]*a_0[i][mu])*theta;
    
    a_0[i][mu]=(r_0*q_0+r_1*q_1+r_2*q_2+r_3*q_3)/k;

    a_1[i][mu]=(-r_0*q_1+r_1*q_0+r_2*q_3-r_3*q_2)/k;

    a_2[i][mu]=(-r_0*q_2-r_1*q_3+r_2*q_0+r_3*q_1)/k;

    a_3[i][mu]=(-r_0*q_3+r_1*q_2-r_2*q_1+r_3*q_0)/k;

    }	

    }

    
}

double heatbath(hb_params_t *hpars)
{
  int isweep;

  beta=hpars->beta;



  for (isweep=0;isweep<hpars->ntraj;isweep++)
    {
      update_heatbath();

 

      /*
        measure();
        if((isweep+1)%(hpars->naccu)==0) {
        print_meas();
        }
      */

    }

  return 0;
}

/*
void heatbath()
{
    int isweep;

    for (isweep=0;isweep<ntraj;isweep++)
    {
     void choose_avector();
    }

}


void choose_avector()
{
    int nx,ny,nz,nt,mu,nu;
    double b,c,r,x,y,theta,phi;
	       	

       for (nx=0;nx<L;nx++)
    {
       for (ny=0;ny<L;ny++)
    {
       for (nz=0;nz<L;nz++)
    {
       for (nt=0;nt<L;nt++)
    {
       for (mu=0;mu<4;mu++)
    {
      
	do {   
               ranlxd(&c,1);

	       ranlxd(&r,1);

               for (nu=0;nu<4;nu++){

            make_k(a_0[(nx+d[0][mu])%L][(ny+d[1][mu])%L][(nz+d[2][mu])%L][(nt+d[3][mu])%L][nu],a_1[(nx+d[0][mu])%L][(ny+d[1][mu])%L][(nz+d[2][mu])%L][(nt+d[3][mu])%L][nu],a_2[(nx+d[0][mu])%L][(ny+d[1][mu])%L][(nz+d[2][mu])%L][(nt+d[3][mu])%L][nu],a_3[(nx+d[0][mu])%L][(ny+d[1][mu])%L][(nz+d[2][mu])%L][(nt+d[3][mu])%L][nu],a_0[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu],a_1[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu],a_2[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu],a_3[(nx+d[0][nu])%L][(ny+d[1][nu])%L][(nz+d[2][nu])%L][(nt+d[3][nu])%L][mu],p_0[nx][ny][nz][nt][mu][nu],p_1[nx][ny][nz][nt][mu][nu],p_2[nx][ny][nz][nt][mu][nu],p_3[nx][ny][nz][nt][mu][nu],a_0[nx][ny][nz][nt][nu],a_1[nx][ny][nz][nt][nu],a_2[nx][ny][nz][nt][nu],a_3[nx][ny][nz][nt][nu],q_0[nx][ny][nz][mu][nu],q_1[nx][ny][nz][nt][mu][nu],q_2[nx][ny][nz][mu][nu],q_3[nx][ny][nz][nt][mu][nu],a_0[(nx+d[0][mu]-d[0][nu])%L][(ny+d[1][mu]-d[0][nu])%L][(nz+d[2][mu]-d[2][nu])%L][(nt+d[3][mu]-d[3][nu])%L][nu],a_1[(nx+d[0][mu]-d[0][nu])%L][(ny+d[1][mu]-d[1][nu])%L][(nz+d[2][mu]-d[2][nu])%L][(nt+d[3][mu]-d[3][nu])%L][nu],a_2[(nx+d[0][mu]-d[0][nu])%L][(ny+d[1][mu]-d[0][nu])%L][(nz+d[2][mu]-d[2][nu])%L][(nt+d[3][mu]-d[3][nu])%L][nu],a_3[(nx+d[0][mu]-d[0][nu])%L][(ny+d[1][mu]-d[1][nu])%L][(nz+d[2][mu]-d[2][nu])%L][(nt+d[3][mu]-d[3][nu])%L][nu],a_0[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],a_1[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],a_2[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],a_3[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],r_0[nx][ny][nz][nt][mu][nu],r_1[nx][ny][nz][nt][mu][nu],r_2[nx][ny][nz][nt][mu][nu],r_3[nx][ny][nz][nt][mu][nu],a_1[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],a_2[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][mu],a_3[(nx-d[0][nu])%L][(ny-d[1][nu])%L][(nz-d[2][nu])%L][(nt-d[3][nu])%L][nu],s_0[nx][ny][nz][nt][mu][nu],s_1[nx][ny][nz][nt][mu][nu],s_2[nx][ny][nz][nt][mu][nu],s_3[nx][ny][nz][nt][mu][nu]);

	       }

	 k_0[nx][ny][nz][nt][mu]=q_0[nx][ny][nz][nt][mu][0]+s_0[nx][ny][nz][nt][mu][0]+q_0[nx][ny][nz][nt][mu][1]+s_0[nx][ny][nz][nt][mu][1]+q_0[nx][ny][nz][nt][mu][2]+s_0[nx][ny][nz][nt][mu][2]+q_0[nx][ny][nz][nt][mu][3]+s_0[nx][ny][nz][nt][mu][3]-q_0[nx][ny][nz][nt][mu][mu]-s_0[nx][ny][nz][nt][mu][mu];

	 k_1[nx][ny][nz][nt][mu]=q_1[nx][ny][nz][nt][mu][0]+s_1[nx][ny][nz][nt][mu][0]+q_1[nx][ny][nz][nt][mu][1]+s_1[nx][ny][nz][nt][mu][1]+q_1[nx][ny][nz][nt][mu][2]+s_1[nx][ny][nz][nt][mu][2]+q_1[nx][ny][nz][nt][mu][3]+s_1[nx][ny][nz][nt][mu][3]-q_1[nx][ny][nz][nt][mu][mu]-s_1[nx][ny][nz][nt][mu][mu];

         k_2[nx][ny][nz][nt][mu]=q_2[nx][ny][nz][nt][mu][0]+s_2[nx][ny][nz][nt][mu][0]+q_2[nx][ny][nz][nt][mu][1]+s_2[nx][ny][nz][nt][mu][1]+q_2[nx][ny][nz][nt][mu][2]+s_2[nx][ny][nz][nt][mu][2]+q_2[nx][ny][nz][nt][mu][3]+s_2[nx][ny][nz][nt][mu][3]-q_2[nx][ny][nz][nt][mu][mu]-s_2[nx][ny][nz][nt][mu][mu];

	 k_3[nx][ny][nz][nt][mu]=q_3[nx][ny][nz][nt][mu][0]+s_3[nx][ny][nz][nt][mu][0]+q_3[nx][ny][nz][nt][mu][1]+s_3[nx][ny][nz][nt][mu][1]+q_3[nx][ny][nz][nt][mu][2]+s_3[nx][ny][nz][nt][mu][2]+q_3[nx][ny][nz][nt][mu][3]+s_3[nx][ny][nz][nt][mu][3]-q_3[nx][ny][nz][nt][mu][mu]-s_3[nx][ny][nz][nt][mu][mu];

	 k[nx][ny][nz][nt][mu]=sqrt(k_0[nx][ny][nz][nt][mu]*k_0[nx][ny][nz][nt][mu]+k_1[nx][ny][nz][nt][mu]*k_1[nx][ny][nz][nt][mu]+k_2[nx][ny][nz][nt][mu]*k_2[nx][ny][nz][nt][mu]+k_3[nx][ny][nz][nt][mu]*k_3[nx][ny][nz][nt][mu]);
	       
               ranlxd(&c,1);

               b=c*(1-exp(-2*beta*k[nx][ny][nz][nt][mu]))+exp(-2*beta*k[nx][ny][nz][nt][mu]);

               a_0[nx][ny][nz][nt][mu]=1+log(b)/(beta*k[nx][ny][nz][nt][mu]);
		       
	       delta=sqrt(1-a_0[nx][ny][nz][nt][mu]*a_0[nx][ny][nz][nt][mu]);

               ranlxd(&r,1);

      	    } while (delta < r);	

        ranlxd(&x,1);
	ranlxd(&y,1);

        theta=acos(2*x-1);
        phi=PI*y-PI/2;

        a_1[nx][ny][nz][nt][mu]=sqrt(1-a_0[nx][ny][nz][nt][mu]*a_0[nx][ny][nz][nt][mu])*cos(theta)*cos(phi);
        a_2[nx][ny][nz][nt][mu]=sqrt(1-a_0[nx][ny][nz][nt][mu]*a_0[nx][ny][nz][nt][mu])*cos(theta)*sin(phi);
        a_3[nx][ny][nz][nt][mu]=sqrt(1-a_0[nx][ny][nz][nt][mu]*a_0[nx][ny][nz][nt][mu])*sin(theta);


    }
    
    }

    }

    }

    }
	
}
    

*/

