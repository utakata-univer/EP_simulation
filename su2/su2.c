#define CONTROL
#include "su2.h" 
#include "string.h"

static hb_params_t hb_params;
static int seed;

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


static int read_input(char *input)
{
    FILE* fp;

    fp=fopen(input,"r");
    if (fp==NULL) {
	fprintf(stderr, "Cannot open input file %s \n",input);
	exit(1);
    }

    get_val(fp, "seed",        "%i" ,&seed   );
    get_val(fp, "beta",        "%lf" ,&hb_params.beta   );
    get_val(fp, "ntraj",        "%i" ,&hb_params.ntraj  );

   
    printf("PARAMETERS\n");
    printf("L              %i\n", L);
    printf("ntraj          %i\n", hb_params.ntraj);
    printf("beta         %f\n", hb_params.beta);
    printf("END PARAMETERS\n");

    return 0;
}

int main(int argc, char* argv[])
{
    double pla;

    if (argc != 2) {
	fprintf(stderr, "Number of arguments not correct\n");
	fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
	exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1]);

    /* initialize random number generator */
    rlxd_init(1,seed);

    hopping(hop);

    /* initialize link */
    initialize();

    heatbath(&hb_params);

    /* compute the plaquette */
    pla=plaquette();

    printf("plaquette = %e\n",pla);

    return 0;
}

