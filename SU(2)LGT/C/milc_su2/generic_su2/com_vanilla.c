/******************  com_vanilla.c **************************************/
/* Communications routines for the SU3 program
/*   MIMD version 3. */

/*
   This file is machine dependent.
   Version for single processor machines, not vectorized.
   This is a trivial version, where there is only one node, every site
   is on node ...

   The lattice is stored as an array of sites, where "site" is a
   structure defined in "LATDEF".

   initialize_machine() does any machine dependent setup at the
   very beginning.

   Some routines implement basic functions in a machine independent
   way.
   mynode() returns node number of this node.
   numnodes() returns number of nodes
   g_sync() provides a synchronization point for all nodes.
   g_floatsum() sums a floating point number over all nodes.
   g_doublesum() sums a double over all nodes.
   g_complexsum() sums a single precision complex number over all nodes.
   g_dcomplexsum() sums a double precision complex number over all nodes.
   g_floatmax() finds maximum of a floating point number over all nodes.
   g_doublemax() finds maximum of a double over all nodes.
   broadcast_float()  broadcasts a single precision number from
	node 0 to all nodes.
   broadcast_double()  broadcasts a double precision number
   broadcast_complex()  broadcasts a single precision complex number
   broadcast_dcomplex()  broadcasts a double precision complex number
   send_integer() sends an integer to one other node
   receive_integer() receives an integer
   dclock() returns a double precision time, with arbitrary zero
   terminate() kills the job on all processors

   make_nn_gathers()  makes all necessary lists for communications with
   nodes containing neighbor sites.

   field_pointer_at_coordinates() returns the address of a field in
   the lattice given its coordinates.
   field_pointer_at_direction() returns the address of a field in the
   lattice at a direction from a given site.
   cleanup_field_pointer() frees the buffers that field_pointer...
   allocated.

   start_gather() starts asynchronous sends and receives required
   to gather neighbors.
   wait_gather()  waits for receives to finish, insuring that the
   data has actually arrived.
   cleanup_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
   restart_gather() repeats the internode communications of a previous
   gather.

   send_parameters() sends a structure of type params to all nodes.
   get_parameters() receives a structure of type params.
   send_field() sends a field to one other node.
   get_links() receives a field from some other node.
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "globaldefs.h"
#include <su2.h>
#include LATDEF
#include <comdefs.h>
#ifndef _IBMR2
#ifndef __osf__
char *malloc();
#endif
#endif

    /* Maintain lists of headers to lists of comlinks */
    /* One list for sites to be received, one for sites to be sent.
		These will often point to the same things */
	    comlink ** neighborlist, **neighborlist_send;
    /* addresses of neighboring sites, NULL if off-node */
    /* neighbor[X][i] is a pointer to the neighbor of site lattice[i] in
       shuffle number X */
    site *** neighbor;
    /* Number of gathers (mappings) that have been set up */
    int n_gathers;


/* Machine initialization */
initialize_machine(argc,argv) int argc; char **argv; {}


/* Set up "comlink" structures needed by gather routines.
   make_lattice() must be called first. */
void make_nn_gathers(){
int i;
int neighbor_coords_special();

    /* Allocate space for lists of pointers to neighbor sites.
       (NULL if neighbor not on this node) */
    neighbor = (site ***)malloc(NDIRS*sizeof(site **));
    n_gathers=0;

    for(i=XUP;i<=TUP;i++)
	make_gather(neighbor_coords_special,&i,WANT_INVERSE,
	    ALLOW_EVEN_ODD,SWITCH_PARITY);

    /* Sort into the order we want for nearest neighbor gathers,
	so you can use XUP, XDOWN, etc. as argument in calling them. */
    sort_eight_special( neighbor );
}

/* sort a list of eight pointers into the order we want for the
  nearest neighbor gathers:  XUP,YUP,ZUP,TUP,TDOWN,ZDOWN,YDOWN,XDOWN */
sort_eight_special( pt ) void **pt; {
void *tt[8];
register int i;
    for(i=0;i<8;i++)tt[i]=pt[i];
    for(i=XUP;i<=TUP;i++){pt[i]=tt[2*i]; pt[OPP_DIR(i)]=tt[2*i+1];}
}

/* utility function for finding coordinates of neighbor */
/* This version for use by make_gather for nearest neighbor gathers */
neighbor_coords_special( x,y,z,t,dirpt,fb, x2p,y2p,z2p,t2p) 
int x,y,z,t,*dirpt,fb;	/* coordinates of site, direction (eg XUP), and
				"forwards/backwards"  */
int *x2p,*y2p,*z2p,*t2p;	/* pointers to coordinates of neighbor */
{
int dir;
    dir = (fb==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
    *x2p = x; *y2p = y; *z2p = z; *t2p = t;
    switch(dir){
	case XUP: *x2p = (x+1)%nx; break;
	case XDOWN: *x2p = (x+nx-1)%nx; break;
	case YUP: *y2p = (y+1)%ny; break;
	case YDOWN: *y2p = (y+ny-1)%ny; break;
	case ZUP: *z2p = (z+1)%nz; break;
	case ZDOWN: *z2p = (z+nz-1)%nz; break;
	case TUP: *t2p = (t+1)%nt; break;
	case TDOWN: *t2p = (t+nt-1)%nt; break;
	default: printf("BOTCH: bad direction\n"); exit(1);
    }
}

#define RECEIVE 0
#define SEND 1
/* add another gather to the list of tables */
int make_gather( func, args, inverse, want_even_odd, parity_conserve )
int (*func)();		/* function which defines sites to gather from */
int *args;		/* list of arguments, to be passed to function */
int inverse;		/* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
int want_even_odd;	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
int parity_conserve;	/* {SAME,SWITCH,SCRAMBLE}_PARITY */
{
register int i,j,k;	/* scratch */
register site *s;	/* scratch */
int dir;		/* direction */
int x,y,z,t;		/* coordinates */

    /* we will have one or two more gathers */
    if( inverse==WANT_INVERSE ) n_gathers += 2;
    else			n_gathers += 1;
    if(n_gathers > MAX_GATHERS){
	if(this_node==0)printf("Too many gathers = %d\n change MAX_GATHERS\n",
			       n_gathers);
	exit(1);
    }
    /* Allocate more space for lists of pointers to neighbor sites.
       (NULL if neighbor not on this node) */
    neighbor = (site ***)realloc(neighbor, n_gathers*sizeof(site **));
    if( inverse==WANT_INVERSE) {
        neighbor[n_gathers-2] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-2]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-2;	/* index of gather we are working on */
    }
    else {
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("NODE %d: no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-1;
    }

    /* Check to see if mapping has advertised parity and inverse properties */
    /* Also check to see if it returns legal values for coordinates */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
	k=(x+y+z+t)%2;

	if( x<0 || y<0 || z<0  || t<0 || x>=nx || y>=ny || z>=nz || t>=nt){
	    printf("DUMMY! Your gather mapping does not stay in lattice\n");
	    printf("It mapped %d %d %d %d to %d %d %d %d\n",
		s->x,s->y,s->z,s->t,x,y,z,t);
	    terminate(1);
	}
	if( ( parity_conserve==SAME_PARITY &&
		((k==0&&s->parity==ODD)||(k==1&&s->parity==EVEN))
	    )
	  ||( parity_conserve==SWITCH_PARITY &&
		((k==0&&s->parity==EVEN)||(k==1&&s->parity==ODD))
	    ) ){
	    printf("DUMMY! Your gather mapping does not obey claimed parity\n");
	    printf("It mapped %d %d %d %d to %d %d %d %d\n",
		s->x,s->y,s->z,s->t,x,y,z,t);
	    terminate(1);
	}
	if( inverse==OWN_INVERSE ){
	    int x2,y2,z2,t2;
            func( x, y, z, t, args, FORWARDS, &x2,&y2,&z2,&t2);
	    if( s->x!=x2 || s->y!=y2 || s->z!=z2 || s->t!=t2 ){
	        printf(
		    "DUMMY! Your gather mapping is not its own inverse\n");
	        printf("It's square mapped %d %d %d %d to %d %d %d %d\n",
		    s->x,s->y,s->z,s->t,x2,y2,z2,t2);
	        terminate(1);
	    }
	}
    }

    /* RECEIVE LISTS: */
    /* Fill in pointers to sites which are on this node, NULL if
	they are off-node */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, FORWARDS, &x,&y,&z,&t);
        j = node_number(x,y,z,t);	/* node for neighbor site */
        /* if neighbor is on node, set up pointer */
        if( j == mynode() ) neighbor[dir][i]= &(lattice[node_index(x,y,z,t)]);
        else		    neighbor[dir][i]= NULL;
    }

    if( inverse != WANT_INVERSE) return(dir);

    /* INVERSE GATHER */
    /* Now, if necessary, make inverse gather */
    /* In most cases, we can use the same lists as the gather, in one
	form or another.  Of course, by the time you get to here
	you know that inverse = WANT_INVERSE */
    dir++;	/* inverse gather has direction one more than original */

    /* Always set up pointers to sites on this node */
    /* scan sites in lattice */
    FORALLSITES(i,s){
        /* find coordinates of neighbor who sends us data */
        func( s->x, s->y, s->z, s->t, args, BACKWARDS, &x,&y,&z,&t);
        j = node_number(x,y,z,t);	/* node for neighbor site */

        /* if neighbor is on node, set up pointer */
        if( j == mynode() ) neighbor[dir][i]= &(lattice[node_index(x,y,z,t)]);
        else 		    neighbor[dir][i]= NULL;
    }

    return(dir-1);
}

/* utility function for finding coordinates of neighbor */
neighbor_coords( x,y,z,t,dir, x2p,y2p,z2p,t2p) 
int x,y,z,t,dir;	/* coordinates of site, and direction (eg XUP) */
int *x2p,*y2p,*z2p,*t2p;	/* pointers to coordinates of neighbor */
{
    *x2p = x; *y2p = y; *z2p = z; *t2p = t;
    switch(dir){
	case XUP: *x2p = (x+1)%nx; break;
	case XDOWN: *x2p = (x+nx-1)%nx; break;
	case YUP: *y2p = (y+1)%ny; break;
	case YDOWN: *y2p = (y+ny-1)%ny; break;
	case ZUP: *z2p = (z+1)%nz; break;
	case ZDOWN: *z2p = (z+nz-1)%nz; break;
	case TUP: *t2p = (t+1)%nt; break;
	case TDOWN: *t2p = (t+nt-1)%nt; break;
	default: printf("BOTCH: bad direction\n"); exit(1);
    }
}

/* Set up interrupt handlers to handle field_pointer routines
   make_lattice() must be called first. */
static msg_request mreqbuf;	/* global so handler can use it too */
void start_handlers(){
void fillfieldrequest();
    /* nothing to do */
}

/* GATHER ROUTINES */
/* start_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_gather() and cleanup_gather() calls.
   In the single processor version, the routine will always return NULL.

   usage:  tag = start_gather( source, size, direction, parity, dest )
   example:
	msg_tag *tag;
	tag = start_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0] );
	  ** do other stuff **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy therof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here.
	  **
	cleanup_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

    Under certain circumstances it is possible to efficiently gather
    a field that has been previously gathered.  This happens when the
    field being gathered has been modified, but the pointers (the
    destination of start_gather() ) have not been modified.  To use
    restart_gather, the original gather must have been waited for but
    not cleaned up.  The usage is:
	msg_tag *tag;
	tag = start_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0] );
	  ** do other stuff, but don't modify tag or gen_pt[0] **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here, but don't modify tag or
	   gen_pt[0].
	   Do modify the source field phi. **
	restart_gather( F_OFFSET(phi), sizeof(su3_vector), XUP,
	    EVEN, gen_pt[0], tag );
	  ** do other stuff **
	wait_gather(tag);
	  ** gen_pt[0][i] now contains the address of the modified phi.
	   The restart-wait may be repeated as often as desired.  **
	cleanup_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

   Internal convention for message types is (schematically):
	( f_offset(src) | sending_node | dir ) + GATHER_BASE
*/

msg_tag * start_gather(field,size,index,parity,dest)
/* arguments */
field_offset field;	/* which field? Some member of structure "site" */
int size;		/* size in bytes of the field (eg sizeof(su3_vector))*/
int index;		/* direction to gather from. eg XUP - index into
			   neighbor tables */
int parity;		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
char ** dest;		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
	case ODD:
	    FORODDSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
    }

    /* return */
    return(NULL);
}

/* Repeat a gather with the same source and destination as a
  previous gather.  The previous gather must have been waited for
  but not cleaned up.  Pointers to sites on the same node are not
  reset, and the same buffers are reused. */
restart_gather(field,size,index,parity,dest,mbuf)
/* arguments */
field_offset field;	/* which field? Some member of structure "site" */
int size;		/* size in bytes of the field (eg sizeof(su3_vector))*/
int index;		/* direction to gather from. eg XUP - index into
			   neighbor tables */
int parity;		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
char ** dest;		/* one of the vectors of pointers */
msg_tag *mbuf;          /* previously returned by start_gather */
{
    /* Nothing to do */
}

void wait_gather(mbuf) msg_tag *mbuf; {
register int i;
    return;
}

void cleanup_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    return;
}

/* GENERAL_GATHER ROUTINES */
/* start_general_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_general_gather() and
   cleanup_general_gather() calls.
   This list contains msg_tags for all receive buffers, followed by
   a msg_tag with msg_id = 0 and msg_buf = NULL, followed by msg_tags
   for all send buffers, followed by a msg_tag with id=0, buf=NULL.
   If no messages at all are required, the routine will return NULL.
   msg_buf=NULL should be a reliable indicator of no message.

   usage:  tag = start_general_gather( source, size, displacement, parity, dest)
   example:
	msg_tag *tag;
	int disp[4]; 
        disp[XUP]=1; disp[YUP]= -1; disp[ZUP] = disp[TUP] = 0;
	tag = start_general_gather( F_OFFSET(phi), sizeof(su3_vector), disp,
	    EVEN, gen_pt[0] );
	  ** do other stuff **
	wait_general_gather(tag);
	  ** gen_pt[0][i] now contains the address of the phi
	   vector (or a copy thereof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it here.
	  **
	cleanup_general_gather(tag);
	  ** subsequent calls will overwrite the gathered fields. but if you
	   don't clean up, you will eventually run out of space **

*/
struct msg_tmp { int node,count; }; /* temporary structure for keeping
	track of messages to be sent or received */
int g_gather_flag=0;	/* flag to tell if general gather in progress */

msg_tag * start_general_gather(field,size,displacement,parity,dest)
/* arguments */
field_offset field;	/* which field? Some member of structure "site" */
int size;		/* size in bytes of the field (eg sizeof(su3_vector))*/
int *displacement;	/* displacement to gather from. four components */
int parity;		/* parity of sites to which we gather.
			   one of EVEN, ODD or EVENANDODD. */
char ** dest;		/* one of the vectors of pointers */
{
/* local variables */
register int i,j,k;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int tsize;		/* size of entry in messages = sizeof(int)+size */
int disp_parity;	/* parity of displacement vector */
int send_parity;	/* parity of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
msg_tag *mbuf;		/* list of message tags, to be returned */
char *tbuf;		/* temporary buffer pointer */
int n_send_msgs, n_recv_msgs;
int type;		/* message type - nread needs its address */


    /* check for gather already in progress */
    if(g_gather_flag!=0){
	fprintf(stderr,"ERROR: node %d, two general_gathers() at once!\n",
	    mynode() );
	exit(1);
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) Make
	list of nodes from whom we expect messages */
    FORSOMEPARITY(i,s,parity){
	if(displacement[XUP]!=0)tx = (s->x + displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y + displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z + displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t + displacement[TUP] + nt)%nt;
			else    tt = s->t;
	dest[i] = F_PT( &lattice[node_index(tx,ty,tz,tt)], field );
    }

    return(NULL);
}

void wait_general_gather(mbuf) msg_tag *mbuf; {
register int i;
    return;
}

void cleanup_general_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    return;
}


/* FIELD POINTER ROUTINES */

/* Return a pointer to a field in the lattice at some coordinates.
   If the site is on this node, just return its address.  If it is
   on another node, make a buffer, get the data, and return address of
   buffer.
   Example:  
	su3_matrix *pt;
	pt = (su3_matrix *)field_pointer_at_coordinates(
	    F_OFFSET(xlink), sizeof(su3_matrix), x,y,z,t );
	... do stuff here ...
	cleanup_field_pointer( (char *)pt );
*/
char * field_pointer_at_coordinates( field, size, x,y,z,t )
/* arguments */
int field;	/* offset of one of the fields in lattice[] */
int size;	/* size of the field in bytes */
int x,y,z,t;	/* coordinates of point to get field from */
{
int index;
    index=node_index(x,y,z,t);
    /* return address */
    return( F_PT( &(lattice[index]), field ) );
}

/* Return a pointer to a field in the lattice at some direction from
   a site on this node. (Usually the "current" site.)
   Works like field_pointer_at_coordinates. */
char * field_pointer_at_direction( field,size, s, direction )
/* arguments */
field_offset field;	/* offset of one of the fields in lattice[] */
int size;	/* size of the field in bytes */
site *s;	/* pointer to a site on this node */
int direction;	/* direction of site's neighbor to get data from.
		   one of XUP,XDOWN,YUP... */
{
    return( F_PT( neighbor[direction][s-lattice], field ) );
}

/* free any buffers allocated by above routines.  Argument is the
address returned by field_pointer...() */
void cleanup_field_pointer(buf) char * buf; {
}

/* SEND AND RECEIVE PARAMETER PACKETS */
send_parameters(buf) params *buf; {
}
get_parameters(buf) params *buf; {
    printf("BOTCH: this never happens\n");
}

/* SEND AND RECEIVE FIELD */
send_field(buf,size,tonode) char *buf; int size,tonode; {
    printf("BOTCH: this never happens\n");
}
get_field(buf,size) params *buf; {
    printf("BOTCH: this never happens\n");
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="Scalar processor";
char * machine_type(){
    return(name);
}

/* Return my node number */
int mynode(){ return(0); }

/* Return number of nodes */
int numnodes(){ return(1); }

/* Synchronize all nodes */
void g_sync(){
}

/* Sum float over all nodes */
void g_floatsum( fpt ) float *fpt; {
}

/* Sum double over all nodes */
void g_doublesum( dpt ) double *dpt; {
}

/* Sum complex over all nodes */
void g_complexsum( cpt ) complex *cpt; {
}

/* Sum double_complex over all nodes */
void g_dcomplexsum( cpt ) double_complex *cpt; {
}

/* Find maximum of float over all nodes */
void g_floatmax( fpt ) float *fpt; {
}

/* Find maximum of double over all nodes */
void g_doublemax( dpt ) double *dpt; {
}

/* Broadcast floating point number from node zero */
void broadcast_float(fpt) float *fpt; {
}

/* Broadcast double precision floating point number from node zero */
void broadcast_double(dpt) double *dpt; {
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(cpt) complex *cpt; {
}

/* Broadcast double precision complex number from node zero */
void broadcast_dcomplex(cpt) double_complex *cpt; {
}

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(tonode,address) int tonode; int *address; {
    printf("BOTCH: this never happens\n");
}

/* Receive an integer from another node */
/* Note we do not check if this was really meant for us */
void receive_integer(address) int *address; {
    printf("BOTCH: this never happens\n");
}

/* Double precision time */
/* This one wraps around after 36 minutes!! It gives the cpu time,
   not the wall clock time */
double dclock(){
long fine;
    fine = clock();
    return( ((double)fine)/1000000.0 );
}

/* version of exit for multinode processes -- kill all nodes */
terminate(status) int status; {
    printf("Termination: status = %d\n",status);
    exit(status);
}
