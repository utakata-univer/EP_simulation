/******************  com_intel.c *****************************************/
/* Communications routines for the SU3 program
/*   MIMD version 3. */

/*
   This file is machine dependent.
   Intel Paragon version

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

   start_general_gather() starts asynchronous sends and receives required
   to gather fields at arbitrary displacement.
   wait_general_gather()  waits for receives to finish, insuring that the
   data has actually arrived, and sets pointers to received data.
   cleanup_general_gather() frees all the buffers that were allocated, WHICH
   MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.

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
char *malloc();

/* Global variables for the communications stuff */
    /* message types for gather encode the direction, the node sending
	the message, and the field offset of the source.  The base is
	added to this so that gather message types are higher than any
	others we use. The direction is in the lowest order bits, then
	the node number, then the field offset. */
    /* For Intel machine, the sending node must be encoded in the
	message type because it is hard to read a message from a
	particular other node, and in a gather nodes often have
	messages coming in from several other nodes. */
    /* for computing message type in gather */
    int mt_nodeshift;	/* number of bits to shift node number */
    int mt_offshift;	/* number of bits to shift field offset */
    int gmt_offshift;	/* number of bits to shift field offset in
	general gather */
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


/* Machine initialization */
#include <sys/types.h>
#include <nx.h>
initialize_machine(argc,argv) int argc; char **argv; {
    flush_to_zero();	/* really turns off all FP traps */
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
	the largest allowed on intel */
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
	neighborlist[n_gathers-1] = NULL;
	neighborlist_send[n_gathers-1] = NULL;
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
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
}
/* The handler for field requests.  It fills the requests, then posts a
   receive for the next request. */
void fillfieldrequest(type,size,node,pid) int type,size,node,pid;{
register char *buf;
register int trapstate;
void fillfieldrequest();

    trapstate=masktrap(1);
    csend( FIELD_REPLY, F_PT( &(lattice[mreqbuf.index]), mreqbuf.field ),
	mreqbuf.size, node, pid);
    hrecv( FIELD_REQUEST, &mreqbuf, sizeof(msg_request), fillfieldrequest );
    masktrap(trapstate);
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
	   vector (or a copy thereof) on the neighbor of site i in the
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
char ** dest;		/* one of the pointer vectors */
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
	mbuf[i].msg_size =  nsites*size;
	mbuf[i].msg_buf = malloc( nsites*size );
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	mbuf[i].msg_id = irecv(
	    GATHER_MSG_TYPE(field,compt->othernode,index),
	    mbuf[i].msg_buf, nsites*size );
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
	mbuf[i].msg_id = isend( GATHER_MSG_TYPE(field,mynode(),index),
	    mbuf[i].msg_buf, nsites*size, compt->othernode, NODE_PID );
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
char ** dest;		/* one of the pointer vectors */
msg_tag *mbuf;		/* previously returned by start_gather */
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
	/* post receive */
	mbuf[i].msg_id = irecv(
	    GATHER_MSG_TYPE(field, mbuf[i].msg_node, index),
	    mbuf[i].msg_buf, mbuf[i].msg_size );
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
	mbuf[i].msg_id = isend( GATHER_MSG_TYPE(field,mynode(),index),
	    mbuf[i].msg_buf, mbuf[i].msg_size, mbuf[i].msg_node, NODE_PID );
    }
}

void wait_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    if(mbuf==NULL)return;
    /* wait for all receive messages */
    for(i=0; mbuf[i].msg_buf != NULL; i++)msgwait( mbuf[i].msg_id );
    /* wait for all send messages */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
	msgwait( mbuf[i].msg_id );
    }
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
int tsize;		/* size of entry in messages = 2*sizeof(int)+size */
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
char ** dest; 		/* one of the vectors of pointers */
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
	    dest[i] = (char *)( &lattice[node_index(tx,ty,tz,tt)]) + field;
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
	if(mbuf[i].msg_buf==NULL){
	    printf("NO ROOM for mbuf, node %d\n",mynode());exit(1);
	}
	/* post receive */
	mbuf[i].msg_id = irecv(
            GENERAL_GATHER_MSG_TYPE(from_nodes[i].node),
            mbuf[i].msg_buf, nsites*tsize );
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
	mbuf[i+n_recv_msgs+1].msg_id = isend(
	    GENERAL_GATHER_MSG_TYPE(this_node),
            mbuf[i+n_recv_msgs+1].msg_buf, nsites*tsize, to_nodes[i].node,
	    NODE_PID );
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
    g_gather_flag=0;
    if(mbuf==NULL)return;
    for(i=0; mbuf[i].msg_buf != NULL; i++){
	msgwait( mbuf[i].msg_id );
	/* set pointers in sites to correct location */
	for(j=0;j<from_nodes[i].count;j++){
	    /* k = index of site on this node, sent in message */
	    k = *(int *)( mbuf[i].msg_buf + j*tsize );
	    tdest[k] = mbuf[i].msg_buf + j*tsize + 2*sizeof(int);
	}
    }
    if( i > 0)free(from_nodes);
    return;
}

void cleanup_general_gather(mbuf) msg_tag *mbuf; {
register int i,i0;
    if(mbuf==NULL)return;
    /* free all receive buffers */
    for(i=0; mbuf[i].msg_buf != NULL; i++)free( mbuf[i].msg_buf );
    /* wait for all send messages, free all send buffers */
    i0=i+1;	/* index of first send message */
    for(i=i0; mbuf[i].msg_buf != NULL; i++){
        msgwait( mbuf[i].msg_id );
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
    node=node_number(x,y,z,t);
    index=node_index(x,y,z,t);
    if( node==mynode() )return( F_PT( &(lattice[index]), field ) );
    /* make buffer */
    buf=malloc(size);
    /* request other node to send data, wait for answer */
    mreq.field=field;
    mreq.size=size;
    mreq.index=index;
    csendrecv( FIELD_REQUEST, &mreq, sizeof(msg_request), node, NODE_PID,
	FIELD_REPLY, buf, size);
    return(buf);
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
    /* if on my node, return address */
    if( neighbor[direction][s-lattice] != NULL ){
	return( F_PT( neighbor[direction][s-lattice], field ) );
    }
    /* otherwise make buffer, find index and other node */
    neighbor_coords(s->x,s->y,s->z,s->t,direction,&x,&y,&z,&t);
    node=node_number(x,y,z,t);
    index=node_index(x,y,z,t);
    buf=malloc(size);
    /* request other node to send data, wait for answer */
    mreq.field=field;
    mreq.size=size;
    mreq.index=index;
    csendrecv( FIELD_REQUEST, &mreq, sizeof(msg_request), node, NODE_PID,
	FIELD_REPLY, buf, size);
    return(buf);
}

/* free any buffers allocated by above routines.  Argument is the
address returned by field_pointer...() */
void cleanup_field_pointer(buf) char * buf; {
    /* if buf is an address in the lattice, leave it alone.  Otherwise
	it was created by malloc and should be freed */
    if( buf < (char *)lattice || buf >= (char *)(lattice+sites_on_node))
	free(buf);
}

/* Copy "size" bytes from "src" to "dest".  This should be assembler
   coded. */
/**
memcpy( dest, src, size ) char *src, *dest; int size; {
register int j;
    for(j=0;j<size;j++)dest[j]=src[j];
}
**/
/**
memcpy( dest, src, size ) char *src, *dest; int size; {
register int j;
    if( j%4==0){ ** transfer words **
	size >>=2;
        for(j=0;j<size;j++)((int *)dest)[j]=((int *)src)[j];
    }
    else{
        for(j=0;j<size;j++)dest[j]=src[j];
    }
}
**/


/* SEND AND RECEIVE PARAMETER PACKETS */
send_parameters(buf) params *buf; {
    csend(PARAM_TYPE,buf,sizeof(params),ALL_NODES,NODE_PID);
}
get_parameters(buf) params *buf; {
    crecv(PARAM_TYPE,buf,sizeof(params));
}

/* SEND AND RECEIVE FIELD */
/* send_field is to be called only by the node doing the sending */
/* get_field is to be called only by the node to which the field was sent */
send_field(buf,size,tonode) char *buf; int size,tonode; {
    csend(FIELD_TYPE,buf,size,tonode,NODE_PID);
}
get_field(buf,size) params *buf; {
    crecv(FIELD_TYPE,buf,size);
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="Paragon";
char * machine_type(){
    return(name);
}

/* Return my node number */
   /* On Intel machine, "mynode()" is already there */

/* Return number of nodes */
   /* On Intel machine, "numnodes()" is already there */

/* Synchronize all nodes */
void g_sync(){
   gsync();
}

/* Sum float over all nodes */
void g_floatsum( fpt ) float *fpt; {
float work;
    gssum(fpt,1,&work);
}

/* Sum double over all nodes */
void g_doublesum( dpt ) double *dpt; {
double work;
    gdsum(dpt,1,&work);
}

/* Sum complex over all nodes */
void g_complexsum( cpt ) complex *cpt; {
/**
float work;
    gssum( &(cpt->real),1,&work);
    gssum( &(cpt->imag),1,&work);
**/
complex work;
    gssum( cpt,2,&work);
}

/* Sum double_complex over all nodes */
void g_dcomplexsum( cpt ) double_complex *cpt; {
double_complex work;
    gdsum( cpt,2,&work);
}

/* Find maximum of float over all nodes */
void g_floatmax( fpt ) float *fpt; {
float work;
    gshigh(fpt,1,&work);
}

/* Find maximum of double over all nodes */
void g_doublemax( dpt ) double *dpt; {
double work;
    gdhigh(dpt,1,&work);
}

/* Broadcast floating point number from node zero */
void broadcast_float(fpt) float *fpt; {
    if(mynode()==0)csend(BROADCAST_FLOAT_TYPE,fpt,sizeof(float),
	ALL_NODES,NODE_PID);
    else crecv(BROADCAST_FLOAT_TYPE,fpt,sizeof(float));
}

/* Broadcast double precision floating point number from node zero */
void broadcast_double(dpt) double *dpt; {
    if(mynode()==0)csend(BROADCAST_DOUBLE_TYPE,dpt,sizeof(double),
	ALL_NODES,NODE_PID);
    else crecv(BROADCAST_DOUBLE_TYPE,dpt,sizeof(double));
}

/* Broadcast single precision complex number from node zero */
void broadcast_complex(cpt) complex *cpt; {
    if(mynode()==0)csend(BROADCAST_COMPLEX_TYPE,cpt,sizeof(complex),
	ALL_NODES,NODE_PID);
    else crecv(BROADCAST_COMPLEX_TYPE,cpt,sizeof(complex));
}

/* Broadcast double precision complex number from node zero */
void broadcast_dcomplex(cpt) double_complex *cpt; {
    if(mynode()==0)csend(BROADCAST_DCOMPLEX_TYPE,cpt,sizeof(double_complex),
	ALL_NODES,NODE_PID);
    else crecv(BROADCAST_DCOMPLEX_TYPE,cpt,sizeof(double_complex));
}

/* Send an integer to one other node */
/* This is to be called only by the node doing the sending */
void send_integer(tonode,address) int tonode; int *address; {
    csend(SEND_INTEGER_TYPE,address,sizeof(int),tonode,NODE_PID);
}

/* Receive an integer from another node */
/* Note we do not check if this was really meant for us */
void receive_integer(address) int *address; {
    crecv(SEND_INTEGER_TYPE,address,sizeof(int));
}

/* Double precision time */
/* dclock() is an Intel library function */

/* version of exit for multinode processes -- kill all nodes */
#define ALL_PROCS -1
terminate(status) int status; {
    printf("Termination: node %d, status = %d\n",this_node,status);
    killcube(ALL_NODES,ALL_PROCS);
    exit(status);
}
