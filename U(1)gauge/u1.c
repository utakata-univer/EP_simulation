

#define CONTROL
#include "phi4.h"
#include "string.h"


/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initialization of ranlux
 */

static hmc_params_t hmc_params;
static act_params_t act_params;
static int seed;


/**********************************************************************
 *    action
 **********************************************************************/

double action(void)
{
    int nx,ny,nz,nt;
    double S;
    double g =act_params g;

    S=0;

    /* loop over all sites */
    for (nx=1;nx<L+1;nx++)
    {
    for (ny=1;ny<L+1;ny++)
    {
    for (nz=1;nz<L+1;nz++)
    {
    for (nt=1;nt<L+1;nt++)
    {	    
	/*sum over neighbors in positive direction*/
    
	S+=(2-2*cos(a[nx][ny][nz][nt][0]+a[nx+1][ny][nz][nt][1]-a[nx][ny+1][nz][nt][0]-a[nx][ny][nz][nt][1])+(2-2*cos(a[nx][ny][nz][nt][0]+a[nx+1][ny][nz][nt][2]-a[nx][ny][nz+1][nt][0]-a[nx][ny][nz][nt][2]))+(2-2*cos(a[nx][ny][nz][nt][0]+a[nx+1][ny][nz][nt][3]-a[nx][ny][nz][nt+1][0]-a[nx][ny][nz][nt][3])+(2-2*cos(a[nx][ny][nz][nt][1]+a[nx][ny+1][nz][nt][2]-a[nx][ny][nz+1][nt][1]-a[nx][ny][nz][nt][2]))+(2-2*cos(a[nx][ny][nz][nt][1]+a[nx][ny+1][nz][nt][3]-a[nx][ny][nz][nt+1][1]-a[nx][ny][nz][nt][3]))+(2-2*cos(a[nx][ny][nz][nt][2]+a[nx][ny][nz+1][nt][3]-a[nx][ny][nz][nt+1][2]-a[nx][ny][nz][nt][3]))/(2*g*g);
    

    }

    }

    }

    }

    return S;
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

    get_val(fp, "g",      "%lf",&act_params.g );
    get_val(fp, "ntherm",      "%i" ,&hmc_params.ntherm );
    get_val(fp, "ntraj",       "%i" ,&hmc_params.ntraj  );
    get_val(fp, "traj_length", "%lf",&hmc_params.tlength);
    get_val(fp, "nstep",       "%i" ,&hmc_params.nstep  );
    get_val(fp, "seed",        "%i" ,&seed   );
    get_val(fp, "naccu",        "%i" ,&hmc_params.naccu   );


   
    printf("PARAMETERS\n");
    printf("L              %i\n", L);
    printf("DIM            %i\n", D);
    printf("g              %f\n", act_params.g);
    printf("ntherm         %i\n", hmc_params.ntherm);
    printf("ntraj          %i\n", hmc_params.ntraj);
    printf("traj_length    %f\n", hmc_params.tlength);
    printf("nstep          %i\n", hmc_params.nstep);
    printf("naccu          %i\n", hmc_params.naccu);
    printf("END PARAMETERS\n");

    return 0;
}
/**********************************************************************
 *     main
 **********************************************************************/

int main(int argc, char* argv[])
{
    double acc;

    if (argc != 2) {
	fprintf(stderr, "Number of arguments not correct\n");
	fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
	exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1]);

    /* initialize random number generator */
    rlxd_init(1,seed);

    /* initialize phi field */
    initialize();

    /*do the updating*/

    acc=hmc(&act_params, &hmc_params);

    printf("ACCRATE %e\n",acc);


    return 0;
}

