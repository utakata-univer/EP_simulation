#define CONTROL
#include "phi42.h"
#include "string.h"


static hmc_params_t hmc_params;
static int seed;


double magnitization(void)
{
    int i,mu;
    double M;
    
    M=0;

    for (i=0;i<V;i++)
    {   
        for (mu=0;mu<D;mu++) M+=phi[hop[i][mu]];
    }

    return M;
}



static int get_val(FILE* fp, char *str, char* fmt, void* val)
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


static int read_input(char *input)
{
    FILE* fp;

    fp=fopen(input,"r");
    if (fp==NULL) {
	fprintf(stderr, "Cannot open input file %s \n",input);
    }

    get_val(fp, "ntherm",       "%i" ,&hmc_params.ntherm );
    get_val(fp, "ntraj",        "%i" ,&hmc_params.ntraj  );
    get_val(fp, "traj_length",  "%lf",&hmc_params.tlength);
    get_val(fp, "nstep",        "%i" ,&hmc_params.nstep  );
    get_val(fp, "seed",         "%i" ,&seed   );



    printf("PARAMETERS\n");
    printf("L              %i\n", L);
    printf("DIM            %i\n", D);
    printf("ntherm         %i\n", hmc_params.ntherm);
    printf("ntraj          %i\n", hmc_params.ntraj);
    printf("traj_length    %f\n", hmc_params.tlength);
    printf("nstep          %i\n", hmc_params.nstep);
    printf("END PARAMETERS\n");

    return 0;
}

int main(int argc, char* argv[])
{
    double mag;

    if (argc !=2) {
        fprintf(stderr, "Number of arguments not correct\n");
	fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
	exit(1);
    }

    read_input(argv[1]);

    hopping(hop);

    mag=magnitization();

    printf("Magnetization = %17.7e\n",mag);

    return 0;
}

