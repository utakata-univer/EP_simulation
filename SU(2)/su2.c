/* 
 *   File su2.c
 *
 *   Contains the main program and a few other routines from which the
 *   code to simulate the su2 gauge theory can be built. Routines for reading 
 *   in the main parameters of the action and the algorithm as  well as 
 *   the computation of the action are provided.
 *
 *   double plaquette(void)
 *       This routine computes the average plaquette <P>, 
 *
 *   P(x,mu,nu)=
 *   (1/2)Tr(U(x,mu)U(x+\hat mu,nu)U(x+\hat nu, mu)^\dagger U(x,nu)^\dagger), 
 *
 *   for the global field U(x,mu)=a0(x,mu) I + a_k(x,mu)  i\sigma_k, 
 *   averaged over all possible ones in the lattice.
 *       One can compute the Wilson gauge action,
 *   
 *    S= \beta sum_{x,mu<nu}(1-P(x,mu,nu)) 
 *
 *   simpley by S=\beta {D(D-1)/2} V (1-<P>) with the parameters,
 *   D, V from lattice.h and beta from act_params.
 *
 *   static int get_val(FILE* fp, char *str, char* fmt,  void* val)
 *      Routine which reads one line from the input file.
 *      Format of the lines is <keyword> <value>.
 *      Checks if the keyword in string str matches,
 *      then gets the value according to the format in fmt
 *
 *   static int read_input(char *input)
 *      Parses the input file (format as specified in get_val)
 *      and prints the parameters onto the screen. Currently
 *      it reads the basic values for the action and also for the 
 *      MC and the seed of the random number generator.
 */ 

#define CONTROL
#include "su2.h"
#include "string.h"


/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in su2.h
 *  seed for initialization of ranlux
 */
static heatbath_params_t heatbath_params;
static act_params_t act_params;
static int seed;


/**********************************************************************
 *    plaquette
 **********************************************************************/

double plaquette(void)
{
  int n,m;
  int mu,nu;	

  double tmp_plaquette;
  double ave_plaquette;
  int number_plaquette;

  double p0;
  /* double m0,m1,m2,m3; */

  double u0,u1,u2,u3;
  double v0,v1,v2,v3;

  number_plaquette = D*(D-1)*V/2;

  ave_plaquette = 0.0; 

  /* loop over all plaquettes */
  for(n=0;n<V;n++){
    for(mu=0;mu<D;mu++){
      for(nu=mu+1;nu<D;nu++){

	m=hop[n][mu]; 
	u0= a0[n][mu]*a0[m][nu] - a1[n][mu]*a1[m][nu] - a2[n][mu]*a2[m][nu] - a3[n][mu]*a3[m][nu];
	u1= a0[n][mu]*a1[m][nu] + a1[n][mu]*a0[m][nu] - a2[n][mu]*a3[m][nu] + a3[n][mu]*a2[m][nu];
	u2= a0[n][mu]*a2[m][nu] + a2[n][mu]*a0[m][nu] - a3[n][mu]*a1[m][nu] + a1[n][mu]*a3[m][nu];
	u3= a0[n][mu]*a3[m][nu] + a3[n][mu]*a0[m][nu] - a1[n][mu]*a2[m][nu] + a2[n][mu]*a1[m][nu];

	m=hop[n][nu]; 
	v0= a0[n][nu]*a0[m][mu] - a1[n][nu]*a1[m][mu] - a2[n][nu]*a2[m][mu] - a3[n][nu]*a3[m][mu];
	v1= a0[n][nu]*a1[m][mu] + a1[n][nu]*a0[m][mu] - a2[n][nu]*a3[m][mu] + a3[n][nu]*a2[m][mu];
	v2= a0[n][nu]*a2[m][mu] + a2[n][nu]*a0[m][mu] - a3[n][nu]*a1[m][mu] + a1[n][nu]*a3[m][mu];
	v3= a0[n][nu]*a3[m][mu] + a3[n][nu]*a0[m][mu] - a1[n][nu]*a2[m][mu] + a2[n][nu]*a1[m][mu];

	p0= u0*v0+u1*v1+u2*v2+u3*v3;

	tmp_plaquette=p0;
	ave_plaquette += tmp_plaquette;

	/* printf("P(%i,%i,%i)=%e %e %e %e\n",n,mu,nu,p0,p1,p2,p3); */

      }
    }
  }

  ave_plaquette /= number_plaquette;
  return ave_plaquette;

}


/**********************************************************************
 *    initialize
 **********************************************************************/
/*
void initialize(void)
{
    int n,mu;
    double x,y,z,b,c,r,r2,phi;
    double Pi;

    Pi=4.* atan(1.);
       
       for (n=0;n<V;n++){
	 for (mu=0;mu<D;mu++){

	   do{
        
        ranlxd(&x,1);
        y=2.0*x-1.0;
	   
        r = sqrt(1-y*y);
        ranlxd(&r2,1);
        
              }while(r < r2);

	   a0[n][mu]=y;

	   ranlxd(&z,1);
	   ranlxd(&b,1);
           phi=2*Pi*z;
	   c=2*b-1;
           
	   a1[n][mu]=sqrt(1-y*y)*sqrt(1-c*c)*cos(phi);
	   a2[n][mu]=sqrt(1-y*y)*sqrt(1-c*c)*sin(phi);
	   a3[n][mu]=sqrt(1-y*y)*c;
	 }
       }	 
        
}       
*/	 
void initialize(void)
{
    int n,mu;

       for (n=0;n<V;n++){
	 for (mu=0;mu<D;mu++){    
	   a0[n][mu]=1.0;
	   a1[n][mu]=0.0;
	   a2[n][mu]=0.0;
	   a3[n][mu]=0.0;
	 }
       }

}	    

/**********************************************************************
 *     get_val
 **********************************************************************/

static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
    char c[128];

    if(1!=fscanf(fp,"%s",c))
    {
	fprintf(stderr,"Error reading input file at %s\n",str);
	exit(1);
    }

    if(strcmp(str,c)!=0)
    {
	fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
	exit(1);
    }

    if(1!=fscanf(fp,fmt,val))
    {
	fprintf(stderr,"Error reading input file at %s\n",str);
	fprintf(stderr,"Cannot read value format %s\n",fmt);
	exit(1);
    }

    return 0;

}


/**********************************************************************
 *     read_input 
 **********************************************************************/

static int read_input(char *input)
{
    FILE* fp;

    fp=fopen(input,"r");
    if (fp==NULL) {
	fprintf(stderr, "Cannot open input file %s \n",input);
	exit(1);
    }

    get_val(fp, "beta",       "%lf",&act_params.beta   );
    get_val(fp, "ntherm",     "%i" ,&heatbath_params.ntherm );
    get_val(fp, "nsweep",     "%i" ,&heatbath_params.nsweep );
    get_val(fp, "naccu",      "%i" ,&heatbath_params.naccu  );
    get_val(fp, "seed",       "%i" ,&seed   );

   
    printf("PARAMETERS\n");
    printf("L              %i\n", L);
    printf("DIM            %i\n", D);
    printf("beta           %f\n", act_params.beta);
    printf("ntherm         %i\n", heatbath_params.ntherm);
    printf("nsweep         %i\n", heatbath_params.nsweep);
    printf("naccu          %i\n", heatbath_params.naccu);
    printf("END PARAMETERS\n");

    return 0;
}
/**********************************************************************
 *     main
 **********************************************************************/

int main(int argc, char* argv[])
{
    /* int n,mu; */
/*   double P,S;*/

    if (argc != 2) {
	fprintf(stderr, "Number of arguments not correct\n");
	fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
	exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1]);

    /* initialize random number generator */
    rlxd_init(1,seed);

    /* initialize the nearest neighbor field */
    hopping(hop);

    /* initialize link field */
    initialize();

    /* check the initialized link field
    for(n=0;n<V;n++){
      for(mu=0;mu<D;mu++){
	printf("U(%i,%i)=%e\n",n,mu,a0[n][mu]);
      }
    }
    */

/*   P=plaquette();
    S=(D*(D-1)/2)*V*(1.0-P);
    printf("\n");
    printf("<P>  : %e\n",P);
    printf("1-<P>: %e\n",1.0-P);
    printf("S    : %e\n",S);


    do the updating 
    printf("\n");*/
    heatbath(&act_params,&heatbath_params);
    /*
    printf("\n");
    printf("ACCRATE %e\n",acc);
    */

/*    P=plaquette();
    S=(D*(D-1)/2)*V*(1.0-P);
    printf("\n");
    printf("<P>  : %e\n",P);
    printf("1-<P>: %e\n",1.0-P);
    printf("S    : %e\n",S);*/

    return 0;
}

