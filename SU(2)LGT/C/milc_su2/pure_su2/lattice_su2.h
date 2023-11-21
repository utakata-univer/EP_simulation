/****************************** lattice.h ********************************/

/* include file for the MIMD QCD program, version 4
   This file defines global scalars and the fields in the lattice. */


#ifndef _LATTICE_H
#define _LATTICE_H


#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif
#define PI 3.14159265358979323846
#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03
#define VERSION_NUMBER 59354
#define SITERAND	/* Use site-based random number generators */
#define NO_GAUGE_FIX 30    /* constants for gauge fixing */
#define COULOMB_GAUGE_FIX 31

/* Options which control the layout ( for 2d-planes layout )

   With the "EVENFIRST" option the even sites are listed contiguously
   in the first part of lattice[], and the odd sites in the last part.
   In this version there must be an equal number of even and odd sites
   in each plane - in other words one of the two shortest dimensions must
   be even.
*/
#define EVENFIRST

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;			/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int warms,trajecs,steps,stepsQ,propinterval;
EXTERN	int saveflag;
EXTERN	float beta;
EXTERN	float epsilon;
EXTERN	char startfile[80],savefile[80];
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN	int total_iters;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;





/* =============================================================*/
/* --------------------The lattice -----------------------------*/
/* =============================================================*/
struct site {
/* --------------------The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ---------------------The physical fields, program dependent */
	/* gauge field */
	su2_matrix link[4];  /* connection */

	/* temporary matrices */
	su2_matrix staple;
	su2_matrix ploop;
};
typedef struct site site;
/*===============================================================*/




/* Vectors for addressing */
/* Generic pointers, for gather routines */
EXTERN char ** gen_pt[8];

/* structure for passing simulation parameters to each node */
typedef struct {
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int stepsQ;	/* number of steps for qhb */
	int propinterval;     /* number of trajectories between measurements */
	int startflag;  /* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;   /* what to do with lattice at end */
	float beta;	/* gauge coupling */
	float epsilon;	/* time step */
	char startfile[80],savefile[80];
}  params;

/* macros for "field offset" and "field pointer", used when fields
  are arguments to subroutines */
/* Usage:  fo = F_OFFSET( field ), where "field" is the name of a field
  in lattice.
     address = F_PT( &site , fo ), where &site is the address of the
  site and fo is a field_offset.  Usually, the result will have to be
  cast to a pointer to the appropriate type. (It is naturally a char *).
*/
typedef int field_offset;
#define F_OFFSET(a) \
  ((field_offset)(((char *)&(lattice[0]. a ))-((char *)&(lattice[0])) ))
#define F_PT( site , fo )  ((char *)( site ) + (fo)) 

/* macros to loop over sites of a given parity.
   Usage:  
	int i;
	site *s;
	FOREVENSITES(i,s){
	    commands, where s is a pointer to the current site and i is
	    the index of the site on the node
	}
*/
#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)
#ifdef EVENFIRST
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<even_sites_on_node;i++,s++)
#define FORODDSITES(i,s) \
    for(i=even_sites_on_node,s= &(lattice[i]);i<sites_on_node;i++,s++)
#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]); \
    i< ( (choice)==EVEN ? odd_sites_on_node : sites_on_node); \
    i++,s++)
/**#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]), \
    last_in_loop = ((choice)==EVEN ? odd_sites_on_node : sites_on_node); \
    i< last_in_loop; \
    i++,s++)**/
#else
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if(s->parity==EVEN)
#define FORODDSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if(s->parity==ODD)
#define FORSOMEPARITY(i,s,choice) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if( (s->parity & (choice)) != 0)
#endif	/* end ifdef EVENFIRST */

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

#endif /* _LATTICE_H */

/* converted to su2 */

