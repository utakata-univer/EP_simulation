/******** setup.c *********/
/* MIMD code version 4 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <globaldefs.h>
#include <su2.h>
#include "lattice.h"
#include <comdefs.h>

/* Each node has a params structure for passing simulation parameters */
params par_buf;

int  setup()   {
   void prnset(),neighset();
   void make_lattice(),make_nn_gathers(),start_handlers(),setup_layout();
   int prompt;

        /* print banner, get volume, nflavors, seed */
    prompt=initial_set();
        /* initialize the node random number generator */
    initialize_prn(&node_prn,iseed,volume+mynode());
        /* Initialize the layout functions, which decide where sites live */
    setup_layout();
        /* allocate space for lattice, set up coordinate fields */
    make_lattice();
        /* set up neighbor pointers and comlink structures */
    make_nn_gathers();
        /* start interrupt handler for field_pointer() requests */
    start_handlers();

    return(prompt);
}


/* SETUP ROUTINES */
int initial_set(){
int prompt,status,get_i();
float get_f();
    /* On node zero, read lattice size, seed, nflavors and send to others */
    if(mynode()==0){
        /* print banner */
        printf("Pure gauge SU2\n");
        printf("Overrelaxed/quasi-heat bath algorithm\n");
        printf(
            "type 0 for no prompts, 1 for prompts or 2 for list of prompts\n");
        status=getprompt(&prompt);
        if (status != 0) 
            {printf("error in input: initial prompt\n");return(-1);}
        nx=get_i(prompt,"nx");
        ny=get_i(prompt,"ny");
        nz=get_i(prompt,"nz");
        nt=get_i(prompt,"nt");
        printf("lattice dimensions = %d %d %d %d\n",nx,ny,nz,nt);
        iseed=get_i(prompt,"iseed");
        printf("random number seed = %d\n",iseed);

        /* fill part of parameter buffer and send it */
        par_buf.nx=nx;
        par_buf.ny=ny;
        par_buf.nz=nz;
        par_buf.nt=nt;
        par_buf.iseed=iseed;
        send_parameters(&par_buf);
    } /* end if(mynode()==0) */
    else {
        get_parameters(&par_buf);
        nx=par_buf.nx;
        ny=par_buf.ny;
        nz=par_buf.nz;
        nt=par_buf.nt;
        iseed=par_buf.iseed;
    }
    
    this_node = mynode();
    number_of_nodes = numnodes();
    volume=nx*ny*nz*nt;
    return(prompt);
}

void make_lattice(){
register int i,j;               /* scratch */
int x,y,z,t;            /* coordinates */
    /* allocate space for lattice, fill in parity, coordinates and index.  */
    lattice = (site *)malloc( sites_on_node * sizeof(site) );
    if(lattice==NULL){
        printf("NODE %d: no room for lattice\n",this_node);
        terminate(1);
    }
   /* Allocate address vectors */
    for(i=0;i<8;i++){
        gen_pt[i] = (char **)malloc(sites_on_node*sizeof(char *) );
        if(gen_pt[i]==NULL){
            printf("NODE %d: no room for pointer vector\n",this_node);
            terminate(1);
        }
    }

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
        if(node_number(x,y,z,t)==mynode()){
            i=node_index(x,y,z,t);
            lattice[i].x=x;     lattice[i].y=y; lattice[i].z=z; lattice[i].t=t;
            lattice[i].index = x+nx*(y+ny*(z+nz*t));
            if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
            else                 lattice[i].parity=ODD;
#ifdef SITERAND
            initialize_prn( &(lattice[i].site_prn) , iseed, lattice[i].index);
#endif
        }
    }
}

/* read in parameters and coupling constants    */
int readin(prompt) int prompt;  {
/* read in parameters for su2 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input  */

int status;
float x,get_f();
int get_i();
char savebuf[20];
void save_ascii(),restore_ascii(),save_binary(),restore_binary();
void coldlat();
double dtime;

    /* On node zero, read parameters and send to all other nodes */
    if(this_node==0){

        printf("\n\n");
    
        /* warms, trajecs */
        warms=get_i(prompt,"warms");
        trajecs=get_i(prompt,"trajecs");
        printf("warms = %d trajecs = %d\n",warms,trajecs);
    
        /* trajectories between propagator measurements */
        propinterval=get_i(prompt,"traj_between_meas");
        printf("no. of traj. between prop. meas. = %d\n",propinterval);
    
        /* get couplings and broadcast to nodes */
        /* beta */
        beta=get_f(prompt,"beta");
        printf("plaq. beta = %f\n",(double)beta);

#ifdef HMC_ALGORITHM
        /* microcanonical time step */
        epsilon=get_f(prompt,"microcanonical_time_step");
        printf("microcanonical time step = %f\n",(double)epsilon);
#endif  
#ifdef RMD_ALGORITHM
        /* microcanonical time step */
        epsilon=get_f(prompt,"microcanonical_time_step");
        printf("microcanonical time step = %f\n",(double)epsilon);
#endif   
        /*microcanonical steps per trajectory */
        steps=get_i(prompt,"steps_per_trajectory");
        printf("number of microcanonical steps per trajectory = %d\n",steps);
    
#ifdef ORA_ALGORITHM
        /*qhb steps per trajectory */
        stepsQ=get_i(prompt,"qhb_steps");
        printf("number of qhb steps per trajectory = %d\n",stepsQ);
#endif   
    
        if (prompt!=0) printf(
            "enter 'continue', 'fresh', 'reload' or reload_binary\n");
        status=scanf("%s",savebuf);
        if(strcmp("fresh",savebuf) == 0 ){
           startflag = FRESH;
        }
        else if(strcmp("reload",savebuf) == 0 ) {
           startflag = RELOAD_ASCII;
            /*read name of file and load it */
            if(prompt!=0)printf("enter name of file containing lattice\n");
            status=scanf("%s",startfile);
            if(status !=1) {
                printf("error in input: file name read\n"); terminate(1);
            }
        }
        else if(strcmp("reload_binary",savebuf) == 0 ) {
           startflag = RELOAD_BINARY;
            /*read name of file and load it */
            if(prompt!=0)printf("enter name of file containing lattice\n");
            status=scanf("%s",startfile);
            if(status !=1) {
                printf("error in input: file name read\n"); terminate(1);
            }
        }
        else if(strcmp("continue",savebuf) == 0 ) {
           startflag = CONTINUE;
            if(this_node==0)printf("continuing with current lattice\n");
        }
        else{
            printf("error in input: lattice_command is invalid\n"); terminate(1);
        }
 
 
 
#ifdef ORA_ALGORITHM
       if (prompt!=0) printf(
           "enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
       status=scanf("%s",savebuf);
       if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
          fixflag = COULOMB_GAUGE_FIX;
          if(this_node==0)printf("fixing to coulomb gauge\n");
       }
       else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
          fixflag = NO_GAUGE_FIX;
          if(this_node==0)printf("NOT fixing the gauge\n");
       }
       else{
           printf("error in input: fixing_command is invalid\n"); terminate(1);
       }
#endif
 
 

        
        if (prompt!=0) printf(
            "'forget' lattice at end,  'save', or 'save_binary'\n");
        status=scanf("%s",savebuf);
        if(status !=1) {
            printf("error in input: 'save' or 'forget'\n");
            terminate(1);
        }
        if(strcmp("save",savebuf) == 0 )  {
            if(prompt!=0)printf("enter filename\n");
            status=scanf("%s",savefile);
            if(status !=1){
                printf("error in input: save filename\n"); terminate(1);
            }
            printf("lattice to be saved in %s\n",savefile);
            saveflag=SAVE_ASCII;
        }
        else if(strcmp("save_binary",savebuf) == 0 ) {
            if(prompt!=0)printf("enter filename\n");
            status=scanf("%s",savefile);
            if(status !=1){
                printf("error in input: save filename\n"); terminate(1);
            }
            printf("lattice to be saved in %s\n",savefile);
            saveflag=SAVE_BINARY;
        }
        else saveflag=FORGET;

        /* make parameter structure and send it */
        par_buf.warms = warms;
        par_buf.trajecs = trajecs;
        par_buf.steps = steps;
        par_buf.stepsQ = stepsQ;
        par_buf.propinterval = propinterval;
        par_buf.startflag = startflag;
        par_buf.fixflag = fixflag;
        par_buf.saveflag = saveflag;
        par_buf.epsilon = epsilon;
        par_buf.beta = beta;
        strcpy(par_buf.startfile,startfile);
        strcpy(par_buf.savefile,savefile);
        send_parameters(&par_buf);
    } /* end if(this_node==0) */
    else{
        get_parameters(&par_buf);
        warms = par_buf.warms;
        trajecs = par_buf.trajecs;
        steps = par_buf.steps;
        stepsQ = par_buf.stepsQ;
        propinterval = par_buf.propinterval;
        startflag = par_buf.startflag;
        fixflag = par_buf.fixflag;
        saveflag = par_buf.saveflag;
        epsilon = par_buf.epsilon;
        beta = par_buf.beta;
        strcpy(startfile,par_buf.startfile);
        strcpy(savefile,par_buf.savefile);
    }
if(beta<0)exit(0);

    dtime = -dclock();
    /* Do whatever is needed to get lattice */
    switch(startflag){
        case CONTINUE:  /* do nothing */
            break;
        case FRESH:     /* cold lattice */
            coldlat();
            break;
        case RELOAD_ASCII:      /* read Ascii lattice */
            restore_ascii(startfile,beta,0.0);
            break;
        case RELOAD_BINARY:     /* read binary lattice */
            restore_binary(startfile,beta,0.0);
            break;
        default:
            if(this_node==0)printf("Bad startflag %d\n",startflag);
            terminate(1);
    }
        dtime += dclock();
        if(this_node==0){
            printf("Time for getting lattice = %e seconds\n",dtime);
        }
    return(0);
}

/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string preceed the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */
/* getprompt gets the initial value of prompt */

float get_f(prompt,variable_name_string)
        int prompt;
        char *variable_name_string;
{
int s;
float x;
char checkname[80];
        if(prompt)  {
                printf("enter %s ",variable_name_string);
                s=scanf("%e",&x);
                if(s == 1) return(x);
        }
        else  {
                s=scanf("%s%e",checkname,&x);
                if (s == EOF) terminate(0);
                if(s == 2 && strcmp(checkname,variable_name_string) == 0)
                   return(x);
        }
        printf("error in input: %s\n",variable_name_string);
        terminate(1);
}
int get_i(prompt,variable_name_string)
        int prompt;
        char *variable_name_string;
{
int s,i;
char checkname[80];
        if(prompt)  {
                printf("enter %s ",variable_name_string);
                s=scanf("%d",&i);
                if (s == 1) return(i);
        }
        else  {
                s=scanf("%s%d",checkname,&i);
                if (s == EOF) terminate(0);
                if(s == 2 && strcmp(checkname,variable_name_string) == 0)
                   return(i);
        }
        printf("error in input: %s\n",variable_name_string);
        terminate(1);
}
int getprompt(prompt) int *prompt; {
char initial_prompt[80];
int stat;
void printprompts();

        scanf("%s",initial_prompt);
        if(strcmp(initial_prompt,"prompt") == 0)  {
           stat=scanf("%d",prompt);
           if (stat != 1) return(1);
        }
        else if(strcmp(initial_prompt,"0") == 0)
           *prompt=0;
        else if(strcmp(initial_prompt,"1") == 0)
           *prompt=1;
        else if(strcmp(initial_prompt,"2") == 0) {
           *prompt=1; printprompts();
        }
        else return(1);

        return(0);
}

void printprompts() {
printf("Here is the list of prompts in the order they are to appear.\n");
printf("A choice among keywords in denoted by {  }.\n");
printf("Optional filenames are only needed if reloading or saving a lattice.\n");
printf("\t prompt\n\t nx\n\t ny\n\t nz\n\t nt\n");
printf("\t iseed\n\t warms\n\t trajecs\n\t traj_between_meas\n");
printf("\t beta\n");
#ifdef RMD_ALGORITHM
printf("\t microcanonical_time_step\n");
#endif
#ifdef HMC_ALGORITHM
printf("\t microcanonical_time_step\n");
#endif
printf("\t steps_per_trajectory\n");
#ifdef ORA_ALGORITHM
printf("\t microcanonical_time_step\n\t qhb_steps\n");
printf("\t'no_gauge_fix', or 'coulomb_gauge_fix'\n");
#endif
printf("\t {continue, fresh, reload, reload_binary}");
printf(" [filename]\n");
printf("\t {forget, save, save_binary}");
printf(" [filename]\n");
}

void coldlat()  {
/* sets link matrices to unit matrices */
  register int i,j,dir;
  register site *sit;

  FORALLSITES(i,sit){
    for(dir=XUP;dir<=TUP;dir++){
      sit->link[dir].e[0] = 1.0;
      for(j=1; j<4; j++)  {
	sit->link[dir].e[j] = 0.0;
      }
    }
  }
  if(this_node==0)printf("new lattice loaded\n");
}
