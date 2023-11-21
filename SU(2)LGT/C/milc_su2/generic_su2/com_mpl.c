/******************  com_mpl.c **************************************/
/* Communications routines for the SU3 program */
/*   MIMD version 4. */

/* C. DeTar port 23 November 1993 */

/*
   This file is machine dependent.
        FOR MPL VERSION ?
        This version implements asynchronous sends and receives.
	No interrupts that I know of, so field_pointer is a problem.

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
   g_floatmax() finds maximum of a floating point number over all nodes.
   g_doublemax() finds maximum of a double over all nodes.
   broadcast_float()  broadcasts a single precision number from
	node 0 to all nodes.
   broadcast_double()  broadcasts a double precision number
   broadcast_complex()  broadcasts a single precision complex number
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
#include <su3.h>
#include LATDEF
#include <comdefs.h>
#include <time.h>   /* For clock() */

#define TRUE 1
#define FALSE 0
#define MAX_WAIT 500000
char mpl_node_name[32] = "Dummy";
int sync_count;         /* diagnostic */
int mpl_tid;		/* task id supplied by MPL  */
int mpl_gid;            /* group identifier supplied by MPL */
int mpl_label;          /* group label (arbitrary integer) */
int mpl_msgid;          /* assigned message number */
int mpl_nbytes;
int mpl_source;
int mpl_type;
int mpl_number_nodes;
int mpl_node_tid[MAX_NUMBER_NODES]; /*  task id number for 
					     logical node number */
int my_logical_node;	/* used as mynode value */

int sum_rcv_wait_count,sum_rcv_wait_count2,rcv_waits;
int sum_snd_wait_count,sum_snd_wait_count2,snd_waits;
float sum_msg_length;
int msgs;

int terminate();

#ifndef _IBMR2
char *malloc();
#endif

/* Global variables for the communications stuff */
    /* message types for gather encode the direction, the node sending
	the message, and the field offset of the source.  The base is
	added to this so that gather message types are higher than any
	others we use. The direction is in the lowest order bits, then
	the node number, then the field offset. */
    /* for computing message type in gather */
    int mt_nodeshift;	/* number of bits to shift node number */
    int mt_offshift;	/* number of bits to shift field offset */
    /* macro to compute the message type */
#define GATHER_MSG_TYPE(offset,node,dir) \
   (GATHER_BASE + ((offset)<<mt_offshift) + ((node)<<mt_nodeshift) + (dir))
#define GENERAL_GATHER_MSG_TYPE(node) (GENERAL_GATHER_BASE + node)

    /* Maintain lists of headers to lists of comlinks */
    /* One list for sites to be received, one for sites to be sent.
	These will often point to the same things */
    comlink ** neighborlist, **neighborlist_send;
    /* addresses of neighboring sites, NULL if off-node */
    /* neighbor[X][i] is a pointer to the neighbor of site lattice[i] in
       gather number X */
    site *** neighbor;
    /* Number of gathers (mappings) that have been set up */
    int n_gathers;

#include <signal.h>
/* Machine initialization */
initialize_machine(argc,argv) int argc; char **argv; {
  int i,j,n,nodes_started,node_tid;
  
  /* Get my task id. */
  if(mpc_environ(&mpl_number_nodes,&mpl_tid))
    {
      fprintf(stderr,"initialize_machine(0): mpc_environ returned %d\n",mperrno);
      mpc_stopall(1); exit(1);
    }
  
  fprintf(stderr,"initialize_machine: node(%d) checking in\n",mpl_tid);
  fflush(stderr);
  
  /* For MPL the task id is the same as the logical node number */
  my_logical_node = mpl_tid;

  if( mpl_tid == 0 ) /* If I am node 0 */
    {
      if(mpl_number_nodes > MAX_NUMBER_NODES) {
	printf("initialize_machine(0): Number of processes exceeds dimension MAX_NUMBER_NODES %d\n",
	       MAX_NUMBER_NODES);
	mpc_stopall(1); exit(1);
      }
      
      /* Once started, control-C should kill all processes */
      signal(SIGHUP,terminate);
      signal(SIGINT,terminate);
      signal(SIGQUIT,terminate);

    } /* End if node 0 */

  /* Set up group consisting of all nodes */
  for(i = 0; i < mpl_number_nodes; i++)
    mpl_node_tid[i] = i;
  mpl_label = 1;
  if(mpc_group(mpl_number_nodes,mpl_node_tid,mpl_label,&mpl_gid))
    {
      fprintf(stderr,"initialize_machine(0): mpc_group returned %d\n",mperrno);
      mpc_stopall(1); exit(1);
    }

  /* Initialize global message statistics */
  sum_rcv_wait_count = 0;
  sum_rcv_wait_count2 = 0;
  rcv_waits = 0;
  sum_snd_wait_count = 0;
  sum_snd_wait_count2 = 0;
  snd_waits = 0;
  sum_msg_length = 0;
  msgs = 0;
}


/* Set up "comlink" structures needed by gather routines.
   make_lattice() must be called first. */
void make_nn_gathers(){
int i;
int neighbor_coords_special();

   /* Set up variables for constructing message types */
   for( mt_nodeshift=0,i=MAX_GATHERS-1; i>0; i>>=1,mt_nodeshift++);
   for( mt_offshift=mt_nodeshift,i=numnodes()-1; i>0; i>>=1,mt_offshift++);
   /* Assume that fields will be ints or bigger, offset divisible by 4 */
    mt_offshift -= 2;
   /* Check for possible overflow of message type.  I think 2^30-1 is
	the largest allowed on ncube */
    /* When this happens, fix the program */
    if( GATHER_MSG_TYPE(sizeof(site),numnodes(),MAX_GATHERS-1) >= 0x3fffffff ){
	if(this_node==0)printf(
	    "Possible overflow of gather message type - Fix the program\n");
	exit(1);
    }

    /* initialize neighborlist[] */
    neighborlist = (comlink **)malloc(NDIRS*sizeof(comlink *));
    neighborlist_send = (comlink **)malloc(NDIRS*sizeof(comlink *));
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
    sort_eight_special( neighborlist );
    sort_eight_special( neighborlist_send );
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
comlink * copy_list_switch();	/* copy linked list of comlinks,
				switching even and odd */
comlink * make_send_receive_list();	/* make linked list of comlinks,
			for lists of sites to be sent or received */
register int i,j,k;	/* scratch */
register site *s;	/* scratch */
int dir;		/* direction */
int x,y,z,t;		/* coordinates */

    /* we will have one or two more gathers */
    if( inverse==WANT_INVERSE ) n_gathers += 2;
    else			n_gathers += 1;
    if(n_gathers > MAX_GATHERS){
	if(this_node==0)printf("Too many gathers! change MAX_GATHERS\n");
	exit(1);
    }
    /* lengthen neighborlist[] */
    neighborlist = (comlink **)realloc(neighborlist,
	n_gathers*sizeof(comlink *));
    neighborlist_send = (comlink **)realloc(neighborlist_send,
	n_gathers*sizeof(comlink *));
    /* Allocate more space for lists of pointers to neighbor sites.
       (NULL if neighbor not on this node) */
    neighbor = (site ***)realloc(neighbor, n_gathers*sizeof(site **));
    if( inverse==WANT_INVERSE) {
	neighborlist[n_gathers-2] = neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-2] = neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-2] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-2]==NULL){
            printf("make_gather(%d): no room for neighbor vector\n",this_node);
            terminate(1);
        }
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("make_gather(%d): no room for neighbor vector\n",this_node);
            terminate(1);
        }
	dir = n_gathers-2;	/* index of gather we are working on */
    }
    else {
	neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-1] = NULL;
        neighbor[n_gathers-1] = (site **)malloc(sites_on_node*sizeof(site *) );
        if(neighbor[n_gathers-1]==NULL){
            printf("make_gather(%d): no room for neighbor vector\n",this_node);
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

    /* make lists of sites which get data from other nodes.  */
    neighborlist[dir] = make_send_receive_list( func, args, want_even_odd,
        FORWARDS, RECEIVE );

    /* SEND LISTS: */
    /* Now make lists of sites to which we send */
    /* Under some conditions, if mapping is its own inverse we can use
	the lists we have already made */
    if( inverse==OWN_INVERSE && want_even_odd==ALLOW_EVEN_ODD
        && parity_conserve==SWITCH_PARITY ){
        /* Under these conditions the send and receive comlinks are the same */
        /* Just plug head of send list into receive list. */
        neighborlist_send[dir]=neighborlist[dir];
    }  /* end code to use same comlinks for send as for receive */
    else if ( inverse==OWN_INVERSE && 
	( want_even_odd==NO_EVEN_ODD ||
	( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SAME_PARITY ) ) ){
        /* Under these conditions, the even list in a send comlink is the
	    odd list in the receive comlink, and vice versa. */
        neighborlist_send[dir] = copy_list_switch( neighborlist[dir] );
    }
    else{
        /* Make new linked list of comlinks for send lists */
        neighborlist_send[dir] = make_send_receive_list( func, args,
	    want_even_odd, FORWARDS, SEND );
    } /* End general case for send lists */

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

    if( want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SWITCH_PARITY ){
        /* Use same comlinks as inverse gather, switching send and receive.
           Nearest neighbor gathers are an example of this case. */
        neighborlist_send[dir]=neighborlist[dir-1];
        neighborlist[dir]=neighborlist_send[dir-1];
    }
    else if( (want_even_odd==ALLOW_EVEN_ODD && parity_conserve==SAME_PARITY)
	|| want_even_odd==NO_EVEN_ODD  ){
        /* make new comlinks, but use same lists as inverse gather, switching
           send and receive, switching even and odd. */
        neighborlist_send[dir] = copy_list_switch( neighborlist[dir-1] );
        neighborlist[dir] = copy_list_switch( neighborlist_send[dir-1] );
    }
    else {  /* general case.  Really only get here if ALLOW_EVEN_ODD
		and SCRAMBLE_PARITY */

        /* RECEIVE LISTS */
        neighborlist[dir] = make_send_receive_list( func, args, want_even_odd,
	    BACKWARDS, RECEIVE );

        /* SEND LISTS: */
        /* Now make lists of sites to which we send */
        neighborlist_send[dir] = make_send_receive_list( func, args,
	    want_even_odd, BACKWARDS, SEND );
    } /* End making new lists for inverse gather */

    return(dir-1);
}


comlink *  make_send_receive_list( func, args, want_even_odd, forw_back,
send_recv )
int (*func)();		/* function which defines sites to gather from */
int *args;		/* list of arguments, to be passed to function */
int want_even_odd;	/* ALLOW_EVEN_ODD or NO_EVEN_ODD */
int forw_back;		/* FORWARDS or BACKWARDS */
int send_recv;		/* SEND or RECEIVE list */
{
register int i,j,k;	/* scratch */
register site *s;	/* scratch */
int x,y,z,t;		/* coordinates */
int parity;		/*if send, parity of site on other node */
			/*if receive, parity of site on this node */
int *ebuf,*obuf;	/* to be malloc'd */
comlink **combuf;	/* to be malloc'd, remember where comlinks are */
register comlink *compt,**comptpt;
comlink *firstpt;

    /* make temporary buffers of numnodes() integers to count numbers of
       even and odd neighbors on each node */
    ebuf = (int *)malloc( numnodes()*sizeof(int) );
    obuf = (int *)malloc( numnodes()*sizeof(int) );
    combuf = (comlink **)malloc( numnodes()*sizeof(comlink *) );

    /* clear neighbor_numbers */
    for(i=0;i<numnodes();i++){
        ebuf[i] = obuf[i] = 0;
    }

    /* scan sites in lattice */
    FORALLSITES(i,s){
        /* find coordinates, node, and parity of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    parity = s->parity;
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    if( (x+y+z+t)%2==0 )parity=EVEN; else parity=ODD;
	}
	j = node_number(x,y,z,t);

        /* if site is off node, increment neighbor_counter */
        if( j != mynode() ){
	    if( send_recv==RECEIVE ){
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD) ebuf[j]++;
	        else   obuf[j]++;
	    }
	    else{
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD) obuf[j]++;
	        else   ebuf[j]++;
	    }
	}
    }

    firstpt=NULL;
    comptpt = &firstpt;
    /* for each neighbor_counter that is nonzero, create a comlink */
    for(j=0;j<numnodes();j++){
        if( j==mynode() )continue;	/* not for local node */
        if( ebuf[j]==0 && obuf[j]==0)continue;
	    /* no neighbors on this node */

        compt = (comlink *)malloc( sizeof(comlink) );
        *comptpt = compt;
        combuf[j] = compt;	/* to make it easy to find again */
        compt->nextcomlink=NULL;	/* currently terminates list */
        compt->othernode=j;
        compt->n_even_connected = ebuf[j];
        compt->n_odd_connected = obuf[j];
        compt->esitelist = (int *)malloc(
	    (ebuf[j]+obuf[j])*sizeof(int) );
        compt->ositelist = (compt->esitelist)+ebuf[j];
	    /*esitelist and ositelist must be filled in later */

        comptpt = &(compt->nextcomlink);	/* linked list, if we
	    extend it this will get address of next comlink. */
    }

    /* clear neighbor_numbers, to be used as counters now */
    for(j=0;j<numnodes();j++){
        ebuf[j]=obuf[j]=0;
    }

    /* scan sites in node again */
    FORALLSITES(i,s){
        /* find coordinates, node, and parity of receiving site */
	if( send_recv==RECEIVE ){
            func( s->x, s->y, s->z, s->t, args, forw_back, &x,&y,&z,&t);
	    parity = s->parity;
	}
	else {  /* SEND */
            func( s->x, s->y, s->z, s->t, args, -forw_back, &x,&y,&z,&t);
	    if( (x+y+z+t)%2==0 )parity=EVEN; else parity=ODD;
	}
	j = node_number(x,y,z,t);

        /* if neighbor is offnode, add to list in appropriate comlink */
        if( j != mynode() ){
	    if(send_recv==RECEIVE ){
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD ){
	            combuf[j]->esitelist[ebuf[j]] = i;
	            ebuf[j]++;
	        }
	        else{
	            combuf[j]->ositelist[obuf[j]] = i;
	            obuf[j]++;
	        }
	    }
	    else { /*SEND*/
	        if( parity==EVEN || want_even_odd==NO_EVEN_ODD ){
	            combuf[j]->ositelist[obuf[j]] = i;
	            obuf[j]++;
	        }
	        else{
	            combuf[j]->esitelist[ebuf[j]] = i;
	            ebuf[j]++;
	        }
	    }
        }
    }
    /* sort the lists of links according to the ordering of their
       even neighbors in the lower numbered node.  The list of sites
       on the lower numbered node is already in order. */
    for(compt=firstpt; compt != NULL; compt=compt->nextcomlink){
        if(compt->othernode > this_node)continue;
	    /* this is lower numbered node, so don't sort */
	if( send_recv==RECEIVE ) i = forw_back;
	else i = -forw_back;
	sort_site_list(compt->n_odd_connected, compt->ositelist,
	    func, args, i);
	sort_site_list(compt->n_even_connected, compt->esitelist,
	    func, args, i);
    }

    /* free temporary storage */
    free(ebuf); free(obuf); free(combuf);
    return(firstpt);
}

comlink * copy_list_switch( old_compt ) comlink *old_compt; {
comlink *firstpt,*compt;
    /* copy a linked list of comlinks, switching even and odd */
    if( old_compt==NULL )return(NULL);
    firstpt = compt = (comlink *)malloc( sizeof(comlink) );
    do{
        compt->othernode=old_compt->othernode;
        compt->n_even_connected = old_compt->n_odd_connected;
        compt->n_odd_connected = old_compt->n_even_connected;
        compt->esitelist = old_compt->ositelist;
        compt->ositelist = old_compt->esitelist;

        if( old_compt->nextcomlink != NULL)
	    compt->nextcomlink = (comlink *)malloc( sizeof(comlink) );
	else compt->nextcomlink = NULL;
	old_compt=old_compt->nextcomlink; 
	compt = compt->nextcomlink;
    } while( old_compt!=NULL );
    return(firstpt);
}

/* sort a list of sites according to the order of the sites on the
   node with which they comunicate */
sort_site_list(n, list, func, args, forw_back)
int n;		/* number of elements in list */
int *list;	/* pointer to list */
int (*func)();	/* function which defines mapping */
int *args;	/* arguments to pass to function */
int forw_back;	/* look forwards or backwards in map */
{
register int j,k,in1,in2,flag;
register site *s;
int x,y,z,t;
    /* bubble sort, if this takes too long fix it later */
    for(j = n-1; j>0; j--){
        flag=0;
	for(k=0; k<j; k++){
	    s = &(lattice[list[k]]);
	    func(s->x,s->y,s->z,s->t,args,forw_back,&x,&y,&z,&t);
	    in1 = node_index(x,y,z,t);
	    s = &(lattice[list[k+1]]);
	    func(s->x,s->y,s->z,s->t,args,forw_back,&x,&y,&z,&t);
	    in2 = node_index(x,y,z,t);
	    if(in1>in2){
		flag=1;
		in1 = list[k];
		list[k]=list[k+1];
		list[k+1]=in1;
	    }
	}
	if(flag==0)break;
    }
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
    /**
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    **/
if(this_node==0)printf("Can't start interrupt handler yet!!\n");
}
/* The handler for field requests.  It fills the requests, then posts a
   receive for the next request. */
void fillfieldrequest(type,size,node,pid) int type,size,node,pid;{
register char *buf;
register int trapstate;
void fillfieldrequest();

if(type != FIELD_REQUEST){
printf("BOTCH node %d from node %d type %d size %d pid %d\n",
mynode(),node,type,size,pid);
printf("BOTCH: request: size = %d, index = %d, offset = %d\n",
mreqbuf.size,mreqbuf.index,mreqbuf.field);
}
    /* Use "masktraps() for simulator, "masktrap()" for real machine */
    /**
    trapstate=masktrap(1);
    csend( FIELD_REPLY, F_PT( &(lattice[mreqbuf.index]), mreqbuf.field ),
	mreqbuf.size, node, pid);
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    masktrap(trapstate);
    **/
    printf("Oops: called fillfieldrequest\n");
}


/* GATHER ROUTINES */
/* start_gather() returns a pointer to a list of msg_tag's, which will
   be used as input to subsequent wait_gather() and cleanup_gather() calls.
   This list contains msg_tags for all receive buffers, followed by
   a msg_tag with msg_id = 0 and msg_buf = NULL, followed by msg_tags
   for all send buffers, followed by a msg_tag with id=0, buf=NULL.
   If no messages at all are required, the routine will return NULL.
   msg_buf=NULL should be a reliable indicator of no message.

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
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
msg_tag *mbuf;		/* list of message tags, to be returned */
char *tbuf;		/* temporary buffer pointer */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    /* figure out how many send and receive messages this gather will
       involve - chase a linked list. */
    for(n_recv_msgs=0, compt= neighborlist[index]; compt != NULL;
	n_recv_msgs++) compt = compt->nextcomlink;
    for(n_send_msgs=0, compt= neighborlist_send[index];
	compt != NULL; n_send_msgs++) compt = compt->nextcomlink;

    /* allocate a buffer for the msg_tags.  This is dynamically allocated
       because there may be an arbitrary number of gathers in progress
       in any direction. */
    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0, compt= neighborlist[index]; compt != NULL;
	i++,compt = compt->nextcomlink){

	/* allocate buffer to receive neighbors */
	switch (parity){
	    case EVEN:
		nsites = compt->n_even_connected;
		sitelist = compt->esitelist;
		break;
	    case ODD:
		nsites = compt->n_odd_connected;
		sitelist = compt->ositelist;
		break;
	    case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
		sitelist = compt->esitelist;
		break;
	}
	mbuf[i].msg_node = compt->othernode;
	mbuf[i].msg_size = nsites*size;
	mbuf[i].msg_buf = malloc( nsites*size );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for mpl the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GATHER_MSG_TYPE(field, mbuf[i].msg_node, index);
	/* Post receive */
	mpl_source = DONTCARE;
	mpl_type = mbuf[i].msg_id;
	mpc_recv((char *)mbuf[i].msg_buf, mbuf[i].msg_size, &mpl_source,
		 &mpl_type, &mbuf[i].mpl_msgid);
	
	/* set pointers in sites to correct location */
	for(j=0;j<nsites;j++){
             dest[sitelist[j]] = mbuf[i].msg_buf + j*size;
	}
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=0;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* set pointers in sites whose neighbors are on this node.  (If all
	neighbors are on this node, this is the only thing done.) */
    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
	    }}
	    break;
	case ODD:
	    FORODDSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
	    }}
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ if(neighbor[index][j] != NULL){
                dest[j] = F_PT(neighbor[index][j],field);
	    }}
	    break;
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Allocate buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
        switch (parity){
            case EVEN:
		nsites = compt->n_odd_connected;
                break;
            case ODD:
		nsites = compt->n_even_connected;
                break;
            case EVENANDODD:
		nsites = compt->n_even_connected + compt->n_odd_connected;
                break;
        }
	tpt=malloc( nsites*size );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i].msg_node=compt->othernode;
	mbuf[i].msg_size=nsites*size;
	mbuf[i].msg_buf=tpt;
	mbuf[i].msg_OK=FALSE;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	  for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	    memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		   size);
	  }
	}
	if( parity==ODD || parity==EVENANDODD ){
	  for(j=0;j<compt->n_even_connected;j++,tpt += size){
	    memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		   size);
	  }
	}
	/* start the send */
	mbuf[i].msg_id = GATHER_MSG_TYPE(field,mynode(),index);
	if(mpc_send( (char *)mbuf[i].msg_buf, mbuf[i].msg_size,
		    mpl_node_tid[mbuf[i].msg_node],
		    mbuf[i].msg_id, &mbuf[i].mpl_msgid))
	  {
	    printf("send_gather(%d): Error %d sending gather to node %d\n",
		   my_logical_node,mperrno,
		   mpl_node_tid[mbuf[i].msg_node]);
	    terminate(1);
	  }
/**	printf("%dSent %d\n",my_logical_node,mbuf[i].msg_size);fflush(stdout);
	printf("(logical %d): Sent gather %d msg %d %d to node %d\n",
	       my_logical_node,i,
	       mbuf[i].msg_id,mpl_msgid,
	       mpl_node_tid[mbuf[i].msg_node]);
	printf("%dSent %d ",my_logical_node,mbuf[i].msg_id); 
	for(j=0;j<6;j++)
	  printf("%f ",((float *)mbuf[i].msg_buf)[j]);
	printf("... (%d bytes)\n",mbuf[i].msg_size); fflush(stdout); **/

      }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }
    /* return */
    return(mbuf);
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
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */
register char *tpt;	/* scratch pointer in buffers */
int nsites;		/* number of sites in this receive or send */
int *sitelist;		/* list of sites in this receive or send */
char *tbuf;		/* temporary buffer pointer */
register comlink *compt;	/* pointer to current comlink */
int n_send_msgs, n_recv_msgs;

    /* figure out how many send and receive messages this gather will
       involve - Use results of previous start_gather(). */
    n_recv_msgs = n_send_msgs = 0;
    if( mbuf != NULL ){
        while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
        while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
      /* use same buffer to receive neighbors */
      mbuf[i].msg_OK = FALSE;
      /* post receive */
      mpl_source = DONTCARE;
      mpl_type = mbuf[i].msg_id;
      mpc_recv((char *)mbuf[i].msg_buf, mbuf[i].msg_size, &mpl_source,
	       &mpl_type, &mbuf[i].mpl_msgid);
    }

    /* for each node whose neighbors I have */
    for(i=n_recv_msgs+1, compt= neighborlist_send[index]; compt != NULL;
	i++,compt = compt->nextcomlink){
	/* Use same buffer to gather data. Remember that when we gather
	   on even sites we must send odd sites, etc. Since the receiving
	   node lists its even sites first, if we are gathering all
	   sites we must do our odd sites first.  */
	tpt=mbuf[i].msg_buf;
	/* gather data into the buffer */
	if( parity==EVEN || parity==EVENANDODD ){
	    for(j=0;j<compt->n_odd_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->ositelist[j]])+field),
		    size);
	    }
	}
	if( parity==ODD || parity==EVENANDODD ){
	    for(j=0;j<compt->n_even_connected;j++,tpt += size){
	        memcpy( tpt, ((char *)(&lattice[compt->esitelist[j]])+field),
		    size);
	    }
	}
	/* start the send */
	mbuf[i].msg_id = GATHER_MSG_TYPE(field,mynode(),index);
	mbuf[i].msg_OK = FALSE;
	if(mpc_send( (char *)mbuf[i].msg_buf, mbuf[i].msg_size,
		    mpl_node_tid[compt->othernode],
		    mbuf[i].msg_id, &mbuf[i].mpl_msgid))
	  {
	    printf("restart_gather(%d): Error %d sending gather to node %d\n",
		   my_logical_node,mperrno,
		   mpl_node_tid[compt->othernode]);
	    terminate(1);
	  }
/**	printf("%dSent %d\n",my_logical_node,mbuf[i].msg_size);fflush(stdout);
	printf("(logical %d): Resent gather %d msg %d %d to node %d\n",
	       my_logical_node,i,
	       mbuf[i].msg_id,mbuf[i].mpl_msgid,
	       mpl_node_tid[compt->othernode]); 
	printf("%dSent %d ",my_logical_node,mbuf[i].msg_id);
	for(j=0;j<6;j++)printf("%f ",
		  ((float *)mbuf[i].msg_buf)[j]);
	printf("... (%d bytes)\n",mbuf[i].msg_size); fflush(stdout);**/
      }
}

void wait_gather(mbuf) msg_tag *mbuf; {
  /* For asynchronous receives with random order of arrival */
  register int i,bufid;
  int n_send_msgs, n_recv_msgs;
  int all_done;
  int wait_count;
  int msg_count;
  int node_pending[MAX_NUMBER_NODES];

  int j; /* debug */
  g_sync();    /* To prevent missequencing of identical message types */
  if(mbuf==NULL)return;

  /* figure out how many send and receive messages this gather
     involved -  */
  n_recv_msgs = n_send_msgs = 0;
  while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
  while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;

  /* wait for all outgoing messages from this node */
  all_done = FALSE;
  wait_count = 0;
  while(!all_done) {
    all_done = TRUE;
    msg_count = 0;
    for(i=n_recv_msgs+1; i<n_recv_msgs+n_send_msgs+1;i++){
      /* If the message hasn't been cleared then check it */
      if(!mbuf[i].msg_OK)
	{
	  if( (mpl_nbytes = mpc_status(mbuf[i].mpl_msgid))
	     == MPL_NOT_COMPLETED )
	    {
	      all_done = FALSE;
	      node_pending[msg_count] = mbuf[i].msg_node;
	      msg_count++;
	    }
	  else
	    /* (We ignore the possibility of an inactive message) */
	    mbuf[i].msg_OK = TRUE;
	}
    }
    wait_count++;
    if(wait_count > MAX_WAIT)
      {
	fprintf(stderr,"wait_gather: send wait time out node %d\n",
	       my_logical_node);fflush(stdout);
	fprintf(stderr,"node %d send wait for nodes",my_logical_node);
	for(i=0; i<msg_count; i++)printf(" %d",node_pending[i]);
	printf("\n");
	/* Try this instead of terminate to see if we get all the output */
	fflush(stdout);fflush(stderr);exit(1);  
      }
  }
  sum_snd_wait_count += wait_count;
  sum_snd_wait_count2 += wait_count*wait_count;
  snd_waits++;
/*  if(wait_count>MAX_WAIT/10){printf("%dsnd_wait %d\n",
			 my_logical_node,wait_count);fflush(stdout);}*/
  
  /* wait for all arriving messages at this node */
  all_done = FALSE;
  wait_count = 0;
  while(!all_done) {
      all_done = TRUE;
      msg_count = 0;
      for(i=0; mbuf[i].msg_buf !=NULL; i++) {
	  /* If we have processed this message already, go on. */
	  if(!mbuf[i].msg_OK) {
	      /* If it has arrived, grab it.  Otherwise, go on. */
	    if( (mpl_nbytes = mpc_status(mbuf[i].mpl_msgid))
	       == MPL_NOT_COMPLETED )
	      {
		all_done = FALSE;
		node_pending[msg_count] = mbuf[i].msg_node;
		msg_count++;
	      }
	    else
	      /* Check off this receive */
	      {
		/* Check byte count for consistency */
		if(mpl_nbytes!=mbuf[i].msg_size)
		  {
		    printf("wait_gather(%d): Got %d bytes, but expected %d in msg %d\n",
		   mpl_nbytes,mbuf[i].msg_size);
		    terminate(1);
		  }
                sum_msg_length += mpl_nbytes;
		msgs++;
		    
/**		printf("(logical %d): Recd gather %d mpl id %d from node %d\n",
		       my_logical_node,i,mbuf[i].mpl_msgid,
		       mpl_source); fflush(stdout);
		printf("%dRecd %d ",my_logical_node,mbuf[i].msg_id);
		for(j=0;j<6;j++)
		  printf("%f ",((float *)mbuf[i].msg_buf)[j]);
		printf("... (%d bytes)\n",mpl_nbytes); fflush(stdout);**/
		
		/* (We ignore the possibility of an inactive message) */
		mbuf[i].msg_OK = TRUE;
	      }
	  }
	}
    wait_count++;
    if(wait_count > MAX_WAIT)
      {
	fprintf(stderr,"wait_gather: recv wait time out node %d\n",
	       my_logical_node);fflush(stdout);
	printf("node %d recv wait for nodes",my_logical_node);
	for(i=0; i<msg_count; i++)printf(" %d",node_pending[i]);
	printf("\n");
	/* Try this instead of terminate to see if we get all the output */
	fflush(stdout);fflush(stderr);exit(1);  
      }
    }
  sum_rcv_wait_count += wait_count;
  sum_rcv_wait_count2 += wait_count*wait_count;
  rcv_waits++;

/*  if(wait_count>MAX_WAIT/10){printf("%drcv_wait %d\n",
			 my_logical_node,wait_count);fflush(stdout);}*/
  return;
}

void cleanup_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /*  free all send buffers */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	free( mbuf[i].msg_buf );
    }
    /* free the msg_tag buffer */
    free(mbuf);
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
struct msg_tmp *to_nodes, *from_nodes;	/* arrays for messages */
int tsize;		/* size of entry in messages = sizeof(int)+size */
char ** tdest;		/* tdest is copy of dest */
/* from_nodes, tsize and tdest are global because they are set in 
   start_general_gather() and used in wait_general_gather().  This
   works because we allow only one general_gather in progress at a
   time. */

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
int disp_parity;	/* parity of displacement vector */
int send_parity;	/* parity of sites that may be sent */
int tx,ty,tz,tt;	/* temporary coordinates */
int othernode;		/* node sent to or received from */
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
    n_recv_msgs = n_send_msgs = 0;
    tsize = 2*sizeof(int)+size;
        /* Use 2*sizeof int so pointer will be aligned to double word */
    tdest = dest;
    /* find parity of sites that may be sent */
    if( (displacement[XUP]+displacement[YUP]+displacement[ZUP]+
	displacement[TUP])%2 == 0)disp_parity=EVEN;
    else disp_parity=ODD;
    switch(parity){
	case EVENANDODD: send_parity=EVENANDODD; break;
	case EVEN:
	    if( disp_parity==EVEN )send_parity=EVEN;
	    else send_parity=ODD;
	    break;
	case ODD:
	    if( disp_parity==EVEN )send_parity=ODD;
	    else send_parity=EVEN;
	    break;
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
	othernode = node_number(tx,ty,tz,tt);
	if( othernode==this_node ){
	    dest[i] = F_PT( &lattice[node_index(tx,ty,tz,tt)], field );
	}
	else{
	    for(j=0;j<n_recv_msgs;j++)if(from_nodes[j].node==othernode)break;
	    if(j < n_recv_msgs){
		from_nodes[j].count++;
	    }
	    else {
	        if(n_recv_msgs==0){
		    from_nodes = (struct msg_tmp *)malloc(
			sizeof(struct msg_tmp) );
		    from_nodes[0].node = othernode;
		    from_nodes[0].count = 1;
		    n_recv_msgs++;
	        }
	        else{
		    from_nodes = (struct msg_tmp *)realloc( from_nodes,
			(n_recv_msgs+1)*sizeof(struct msg_tmp) );
		    from_nodes[j].node = othernode;
		    from_nodes[j].count = 1;
		    n_recv_msgs++;
	        }
	    }
	}
    }

    /* scan sites of parity we are sending, make list of nodes to which
	we must send messages and the number of messages to each. */
    FORSOMEPARITY(i,s,send_parity){
	if(displacement[XUP]!=0)tx = (s->x - displacement[XUP] + nx)%nx;
			else    tx = s->x;
	if(displacement[YUP]!=0)ty = (s->y - displacement[YUP] + ny)%ny;
			else    ty = s->y;
	if(displacement[ZUP]!=0)tz = (s->z - displacement[ZUP] + nz)%nz;
			else    tz = s->z;
	if(displacement[TUP]!=0)tt = (s->t - displacement[TUP] + nt)%nt;
			else    tt = s->t;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    if(j < n_send_msgs){
		to_nodes[j].count++;
	    }
	    else {
	        if(n_send_msgs==0){
		    to_nodes = (struct msg_tmp *)malloc(sizeof(struct msg_tmp));
		    to_nodes[0].node = othernode;
		    to_nodes[0].count = 1;
		    n_send_msgs++;
	        }
	        else{
		    to_nodes = (struct msg_tmp *)realloc( to_nodes,
			(n_send_msgs+1)*sizeof(struct msg_tmp) );
		    to_nodes[j].node = othernode;
		    to_nodes[j].count = 1;
		    n_send_msgs++;
	        }
	    }
	}
    }

    if( n_recv_msgs==0 && n_send_msgs==0)mbuf=NULL;
    else {
	mbuf = (msg_tag *)malloc(
	    (n_recv_msgs+n_send_msgs+2)*sizeof(msg_tag) );
	if(mbuf==NULL){printf("NO ROOM for mbuf, node %d\n",mynode()); exit(1);}
    }

    /* for each node which has neighbors of my sites */
    for(i=0; i<n_recv_msgs; i++){
	/* allocate buffer to receive neighbors */
	nsites = from_nodes[i].count;
	mbuf[i].msg_node = from_nodes[i].node;
	mbuf[i].msg_size = nsites*tsize;
	mbuf[i].msg_buf = malloc( nsites*tsize );
	mbuf[i].msg_OK = FALSE;
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* for mpl the GATHER_MSG_TYPE identifies the message */
	mbuf[i].msg_id = GENERAL_GATHER_MSG_TYPE(from_nodes[i].node);
	/* Post receive */
	mpl_source = DONTCARE;
	mpl_type = mbuf[i].msg_id;
	mpc_recv((char *)mbuf[i].msg_buf, mbuf[i].msg_size, &mpl_source,
		 &mpl_type, &mbuf[i].mpl_msgid);
    }
    /* terminate receive portion of list */
    if(mbuf != NULL){
        mbuf[n_recv_msgs].msg_id=0;
        mbuf[n_recv_msgs].msg_buf=NULL;
    }

    /* for each node whose neighbors I have */
    for(i=0; i<n_send_msgs; i++){
	/* Allocate buffer to gather data. */
	tpt=malloc( to_nodes[i].count*tsize );
        if(tpt==NULL){printf("NO ROOM for tpt, node %d\n",mynode());exit(1);}
	mbuf[i+n_recv_msgs+1].msg_node=to_nodes[i].node;
	mbuf[i+n_recv_msgs+1].msg_size=to_nodes[i].count*tsize;
	mbuf[i+n_recv_msgs+1].msg_buf=tpt;
	mbuf[i+n_recv_msgs+1].msg_OK=FALSE;
    }

    /* reset to_node counters */
    for(i=0;i<n_send_msgs;i++)to_nodes[i].count=0;
    /* gather data into the buffers. Each entry in the buffers consists
	of the index of the site to which the data is sent, followed by
	the actual data */
    FORSOMEPARITY(i,s,send_parity){
	tx = (s->x - displacement[XUP] + nx)%nx;
	ty = (s->y - displacement[YUP] + ny)%ny;
	tz = (s->z - displacement[ZUP] + nz)%nz;
	tt = (s->t - displacement[TUP] + nt)%nt;
	othernode = node_number(tx,ty,tz,tt);
	if( othernode != this_node ){
	    for(j=0;j<n_send_msgs;j++)if(to_nodes[j].node==othernode)break;
	    tpt = mbuf[j+n_recv_msgs+1].msg_buf +
		to_nodes[j].count*tsize;
	    *(int *)tpt = node_index(tx,ty,tz,tt);
	        /* index of site on other node */
	    memcpy( tpt+2*sizeof(int), F_PT(s,field), size);
	    to_nodes[j].count++;
	}
    }

    /* start the sends */
    for(i=0;i<n_send_msgs;i++){
	nsites = to_nodes[i].count;
	mbuf[i+n_recv_msgs+1].msg_id = GENERAL_GATHER_MSG_TYPE(this_node);
	if(mpc_send( (char *)mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize,
		    mpl_node_tid[to_nodes[i].node],
		    mbuf[i+n_recv_msgs+1].msg_id, &mbuf[i+n_recv_msgs+1].mpl_msgid))
	  {
	    printf("start_general_gather(%d): Error %d sending gen gather to node %d\n",
		   my_logical_node,mperrno,
		   mpl_node_tid[to_nodes[i].node]);
	    terminate(1);
	  }
/**	printf("%dSent %d\n",my_logical_node,mbuf[i].msg_size);fflush(stdout);
	printf("(logical %d): Sent gen gather %d msg %d to node %d\n",
	       my_logical_node,i,mbuf[i+n_recv_msgs+1].mpl_msgid,
	       mpl_node_tid[to_nodes[i].node]);
	printf("Message sent ...\n");
	for(j=0;j<nsites*tsize;j++)printf("%x",
		  ((char *)mbuf[i+n_recv_msgs+1].msg_buf)[j]);
	printf("\n"); **/
      }
    /* terminate list */
    if(mbuf != NULL){
        mbuf[n_send_msgs+n_recv_msgs+1].msg_id=0;
        mbuf[n_send_msgs+n_recv_msgs+1].msg_buf=NULL;
    }

    /* free temporary arrays */
    if( n_send_msgs > 0)free(to_nodes);
    /* mark gather in progress and return */
    g_gather_flag=1;
    return(mbuf);
}

void wait_general_gather(mbuf) msg_tag *mbuf; {
  register int i,j,k;
  register comlink *compt;	/* pointer to current comlink */
  int n_send_msgs, n_recv_msgs;
  int all_done;
  int wait_count;

  g_sync();    /* To prevent missequencing of identical message types */
  g_gather_flag=0;
  if(mbuf==NULL)return;

  /* figure out how many send and receive messages this gather
     involved -  */
  n_recv_msgs = n_send_msgs = 0;
  while( mbuf[n_recv_msgs].msg_buf != NULL) n_recv_msgs++;
  while( mbuf[n_recv_msgs+n_send_msgs+1].msg_buf != NULL ) n_send_msgs++;

  /* wait for all outgoing messages from this node */
  all_done = FALSE;
  wait_count = 0;
  while(!all_done) {
    all_done = TRUE;
    for(i=n_recv_msgs+1; i<n_recv_msgs+n_send_msgs+1;i++){
      /* If the message hasn't been cleared then check it */
      if(!mbuf[i].msg_OK)
	{
	  if( (mpl_nbytes = mpc_status(mbuf[i].mpl_msgid))
	     == MPL_NOT_COMPLETED )
	    all_done = FALSE;
	  else
	    /* (We ignore the possibility of an inactive message) */
	    mbuf[i].msg_OK = TRUE;
	}
    }
    wait_count++;
    if(wait_count > MAX_WAIT)
      {
	fprintf(stderr,"wait_general_gather: send wait time out node %d\n",
	       my_logical_node);fflush(stdout);
	terminate(1);
      }
  }

  sum_snd_wait_count += wait_count;
  sum_snd_wait_count2 += wait_count*wait_count;
  snd_waits++;

/*  if(wait_count>MAX_WAIT/10){printf("%dsnd_wait %d\n",
			 my_logical_node,wait_count);fflush(stdout);}*/

  /* for each node which has neighbors of my sites */
  /* wait for all receive messages */
  all_done = FALSE;
  wait_count = 0;
  while(!all_done) {
    all_done = TRUE;
    for(i=0; mbuf[i].msg_buf !=NULL; i++) {
      /* If we have processed this message already, go on. */
      if(!mbuf[i].msg_OK) {
	/* If it has arrived, grab it.  Otherwise, go on. */
	if( mpc_status(mbuf[i].mpl_msgid) == MPL_NOT_COMPLETED )
	  /* This message not received yet, so we are not done */
	  all_done = FALSE;
	else
	  {
	    /* Check off this receive */
	    /* (We ignore the possibility of an inactive message) */
	    mbuf[i].msg_OK = TRUE;
	    /**		  printf("(logical %d): Got gen gather %d msg %d\n",
	      my_logical_node,i,mbuf[i].mpl_msgid);
	      printf("Message received ...\n");
	      for(j=0;j<mbuf[i].msg_size;j++)printf("%x",
	      ((char *)mbuf[i].msg_buf)[j]);
	      printf("\n"); **/
	    /* set pointers in sites to correct location */
	    for(j=0;j<from_nodes[i].count;j++){
	      /* k = index of site on this node, sent in message */
	      k = *(int *)( mbuf[i].msg_buf + j*tsize );
	      tdest[k] = mbuf[i].msg_buf + j*tsize + 2*sizeof(int);
	    }
	  }
      }
    }
    wait_count++;
    if(wait_count > MAX_WAIT)
      {
	fprintf(stderr,"wait_general_gather: recv wait time out node %d\n",
	       my_logical_node);fflush(stdout);
	terminate(1);
      }
  }
  if( i > 0)free(from_nodes);

  sum_rcv_wait_count += wait_count;
  sum_rcv_wait_count2 += wait_count*wait_count;
  rcv_waits++;

  /* if(wait_count>MAX_WAIT/10){printf("%drcv_wait %d\n",
			 my_logical_node,wait_count);fflush(stdout);}*/
  return;
}

void cleanup_general_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /* wait for all send messages, free all send buffers */
    /* In simulator, messages were synchronous so don't wait */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	free( mbuf[i].msg_buf );
    }
    /* free the msg_tag buffer */
    free(mbuf);
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
register int node,index;
register char *buf;
msg_request mreq;
    /* if on my node, return address */
printf("Oops: tried to use field_pointer()\n"); terminate(1);
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
register int node,index,trapstate,id;
int x,y,z,t;
register char *buf;
msg_request mreq;
printf("Oops: tried to use field_pointer()\n"); terminate(1);
}

/* free any buffers allocated by above routines.  Argument is the
address returned by field_pointer...() */
void cleanup_field_pointer(buf) char * buf; {
    /* if buf is an address in the lattice, leave it alone.  Otherwise
	it was created by malloc and should be freed */
    if( buf < (char *)lattice || buf >= (char *)(lattice+sites_on_node))
	free(buf);
}

/* SEND AND RECEIVE PARAMETER PACKETS TO ALL NODES IN THIS PARTITION*/
send_parameters(buf) params *buf; {
  int i;

  for(i=0; i < mpl_number_nodes; i++) if(i != mynode())
    {
      if(mpc_bsend((char *)buf,sizeof(params), 
		  mpl_node_tid[i],PARAM_TYPE))
	{
	  printf("send_parameters(%d): Error %d sending params to node %d\n",
		 my_logical_node,mperrno,mpl_node_tid[i]);
	  terminate(1);
	}
/**      printf("%dSent %d\n",my_logical_node,sizeof(params));fflush(stdout);**/
    }
}
get_parameters(buf) params *buf; {

    mpl_source = DONTCARE;
    mpl_type = PARAM_TYPE;
    if(mpc_brecv((char *)buf,sizeof(params), 
		 &mpl_source, &mpl_type, &mpl_nbytes))
      {
	printf("get_parameters(%d): Error %d receiving params from node %d\n",
	       my_logical_node,mperrno,mpl_source);
	terminate(1);
      }
}

/* SEND AND RECEIVE FIELD */
send_field(buf,size,tonode) char *buf; int size,tonode; {
  
  if(mpc_bsend((char *)buf,size,
	      mpl_node_tid[tonode],FIELD_TYPE))
    {
      printf("send_field(%d): Error %d sending field to node %d\n",
	     my_logical_node,mperrno,mpl_node_tid[tonode]);
      terminate(1);
    }
/**	printf("%dSent %d to %d\n",my_logical_node,size,mpl_node_tid[tonode]);fflush(stdout);**/
}
get_field(buf,size) params *buf;int size; {
  
  mpl_source = DONTCARE;
  mpl_type = FIELD_TYPE;
  if(mpc_brecv((char *)buf, size,
	       &mpl_source, &mpl_type, &mpl_nbytes))
    {
      printf("get_field(%d): Error %d receiving field from node %d\n",
	     my_logical_node,mperrno,mpl_source);
      terminate(1);
    }
/**	printf("%dRecd %d from %d\n",my_logical_node,size,mpl_source);fflush(stdout);**/
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="mpl";
char * machine_type(){
    return(name);
}

/* Return my node number */
mynode()
{
  return my_logical_node;
}

/* Return number of nodes */
numnodes(){
  return(mpl_number_nodes );
}

/* Synchronize all nodes */
void g_sync(){
/**  printf("%dSync %d\n",my_logical_node,sync_count++);fflush(stdout);**/
  mpc_sync(mpl_gid);
}

/* Sum float over all nodes */
void g_floatsum( fpt ) float *fpt; {
  float sum;
  int othernode,type;
  
/*THIS IS INDETERMINATE AT THE MOMENT BECAUSE THINGS ARE ADDED
IN UNKNOWN ORDER */

  mpc_combine(fpt,&sum,sizeof(float),s_vadd,mpl_gid);

  /* Answer goes in fpt */
  *fpt = sum;
}

/* Sum double over all nodes */
void g_doublesum( dpt ) double *dpt; {
  double sum;
  mpc_combine(dpt,&sum,sizeof(double),d_vadd,mpl_gid);

  /* Answer goes in dpt */
  *dpt = sum;
}

/* Sum complex over all nodes */
void c_vadd (complex *in1, complex *in2, complex *out, int *len)
{
  int i,n;
  n = *len/sizeof(complex);
  for (i=0; i<n; i++)
    CADD(in1[i],in2[i],out[i]);
}

void g_complexsum( cpt ) complex *cpt; {
  complex sum;
  mpc_combine(cpt,&sum,sizeof(complex),c_vadd,mpl_gid);

  /* Answer goes in cpt */
  *cpt = sum;
}

/* Find maximum of float over all nodes */
void g_floatmax( fpt ) float *fpt; {
  float high;

  mpc_combine(fpt,&high,sizeof(float),s_vmax,mpl_gid);

  /* Answer goes in fpt */
  *fpt = high;
}


/* Find maximum of double over all nodes */
void g_doublemax( dpt ) double *dpt; {
  double high;
  
  mpc_combine(dpt,&high,sizeof(double),d_vmax,mpl_gid);

  /* Answer goes in dpt */
  *dpt = high;
}


/* Broadcast floating point number from node zero */
void broadcast_float(fpt) float *fpt; {

  mpl_source = 0;
  mpc_bcast(fpt, sizeof(float), mpl_source, mpl_gid);
/**	printf("%dBcst %d\n",my_logical_node,sizeof(float));fflush(stdout);**/
}

/* Broadcast double precision floating point number from node zero */
void broadcast_double(dpt) double *dpt; {

  mpl_source = 0;
  mpc_bcast(dpt, sizeof(double), mpl_source, mpl_gid);
/**	printf("%dBcst %d\n",my_logical_node,sizeof(double));fflush(stdout);**/
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(cpt) complex *cpt; {

  mpl_source = 0;
  mpc_bcast(cpt, sizeof(complex), mpl_source, mpl_gid);
/**	printf("%dBcst %d\n",my_logical_node,sizeof(complex));fflush(stdout);**/

}

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(tonode,address) int tonode; int *address; {

  if(mpc_bsend((int *)address, sizeof(int),
	      mpl_node_tid[tonode], SEND_INTEGER_TYPE))
    {
      printf("send_integer(%d): Error %d sending integer to node %d\n",
	     my_logical_node,mperrno,mpl_node_tid[tonode]);
      terminate(1);
    }
/**	printf("%dSent %d\n",my_logical_node,sizeof(int));fflush(stdout);**/
}

/* Receive an integer from another node */
/* Note we do not check if this was really meant for us */
void receive_integer(address) int *address; {
  mpl_source = DONTCARE;
  mpl_type = FIELD_TYPE;
  if(mpc_brecv((int *)address,sizeof(int),
	       &mpl_source, &mpl_type, &mpl_nbytes))
    {
      printf("receive_integer(%d): Error %d receiving integer from node %d\n",
	     my_logical_node,mperrno,mpl_source);
      terminate(1);
    }
}

/* Double precision time */
/* This one wraps around after 36 minutes!! It gives the cpu time,
   not the wall clock time */
double dclock(){
long fine;
    fine = clock();
    return( ((double)fine)/CLOCKS_PER_SEC);
}

/* version of exit for multinode processes -- kill all nodes */
terminate(status) int status; {
  printf("Termination: node %d, status = %d\n",my_logical_node,status);
  printf("Message send waiting: node %d %d waits, %e +/- %e\n",
	 my_logical_node,snd_waits,((float)sum_snd_wait_count)/snd_waits,
	 sqrt(((double)sum_snd_wait_count2)/snd_waits - 
	      (((double)sum_snd_wait_count)/snd_waits*((double)sum_snd_wait_count)/snd_waits)));
  printf("Message recv waiting: node %d %d waits, %e +/- %e\n",
	 my_logical_node,rcv_waits,((float)sum_rcv_wait_count)/rcv_waits,
	 sqrt(((double)sum_rcv_wait_count2)/rcv_waits - 
	      (((double)sum_rcv_wait_count)/rcv_waits*((double)sum_rcv_wait_count)/rcv_waits)));
  printf("Avg message length %f\n",sum_msg_length/msgs);
  fflush(stdout);fflush(stderr);
  sleep(1); /* Pause for output flushing */
  mpc_stopall(status);
  exit(status);
}

