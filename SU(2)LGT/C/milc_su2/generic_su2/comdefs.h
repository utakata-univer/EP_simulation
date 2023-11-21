/************************* comdefs.h *************************************/

/* Definitions for communications for the SU3 program on the Intel machine,
   version 4.

   Communications routines will assume that the lattice is stored as an
   array of structures of type "site".
*/

/* message types  (not all are used on all machines) */
#define PARAM_TYPE 11	/* type of parameter message to nodes */
#define FIELD_TYPE 12	/* type of field sent from one node to another */
#define BROADCAST_FLOAT_TYPE 13	/* broadcast of floating point number */
#define BROADCAST_DOUBLE_TYPE 14 /* broadcast of double */
#define BROADCAST_COMPLEX_TYPE 15 /* broadcast of single precision complex */
#define BROADCAST_DCOMPLEX_TYPE 16 /* broadcast of double precision complex */
#define SEND_INTEGER_TYPE 17	/* send an integer to one other node */
#define CM_GATHER_TYPE 18	/* type for CM5 (cooperative) messages */
#define SYNC_TYPE 50		/* Synchronize all nodes */
#define SUM_FLOAT_TYPE 51	/* Sum float over all nodes */
#define SUM_DOUBLE_TYPE 52	/* Sum double over all nodes */
#define SUM_COMPLEX_TYPE 53	/* Sum complex over all nodes */
#define SUM_DCOMPLEX_TYPE 54	/* Sum double_complex over all nodes */
#define MAX_FLOAT_TYPE 55	/* Maximum float over all nodes */
#define MAX_DOUBLE_TYPE 56	/* Maximum double over all nodes */

/* Added for pvm */
#ifdef PVM
#define ANY_MSG -1              /* Any message */

/* Message structures for communication between host and nodes */
/* The first two fields must always be the same for each type */
/* The basic structure must be the shortest */

/* For most messages */
struct hcs_basic {
  int msg_type;
  int node;	/* integer identifies caller's instance */
  int arg1,arg2,arg3;  /* Use depends on which routine */
} ;
#define HCS_BASIC_SIZE (sizeof(struct hcs_basic))
#ifdef PROTO
int put_hcs_basic( struct hcs_basic * hcs );
int get_hcs_basic( struct hcs_basic * hcs );
#endif

/* For printf and scanf calls */
#define STRING_LENGTH 256
struct hcs_stdio {
  int msg_type;
  int node;	/* integer identifies caller's instance */
  int length;
  char s[STRING_LENGTH];
} ;
#define HCS_STDIO_SIZE (sizeof(struct hcs_stdio))
#ifdef PROTO
int put_hcs_stdio( struct hcs_stdio * hcs );
int get_hcs_stdio( struct hcs_stdio * hcs );
#endif

/* For initialization call */
#define MAX_NUMBER_NODES 8
#define HOST_NAME_LENGTH 128
struct hcs_ident {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_instance[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident))
#ifdef PROTO
int put_hcs_ident( struct hcs_ident * hcs );
int get_hcs_ident( struct hcs_ident * hcs );
#endif

union {
  struct hcs_basic basic;
  struct hcs_stdio stdio;
  struct hcs_ident ident;
} hcs;

#define HOST_CALL 77	/* pvm message type for call to host for service */
#define HOST_REPLY 87

/* Message subtypes internal to host-node service calls */
/* (Not used by pvm to identify messages) */
#define PRINTF_HOST_CALL 11
#define SCANF_HOST_CALL 12
#define FPRINTF_HOST_CALL 13 /* Not used */ 
#define FSCANF_HOST_CALL 14 /* Not used */
#define FFLUSH_HOST_CALL 20
#define FOPEN_HOST_CALL 30 /* Not used */
#define FCLOSE_HOST_CALL 31 /* Not used */
#define OPEN_HOST_CALL 32 /* Not used */
#define CLOSE_HOST_CALL 33 /* Not used */
#define CREAT_HOST_CALL 34 /* Not used */
#define READ_HOST_CALL 40 /* Not used */
#define WRITE_HOST_CALL 41 /* Not used */
#define NODE_IDENT_CALL 78	
#define NODES_DONE_HOST_CALL 99

/* end of pvm additions */
#endif	/* end ifdef PVM */

/* Added for pvm */
#ifdef PVM24
#define ANY_MSG -1              /* Any message */

/* Message structures for communication between host and nodes */
/* The first two fields must always be the same for each type */
/* The basic structure must be the shortest */

/* For initialization call */
#define MAX_NUMBER_NODES 8
#define HOST_NAME_LENGTH 128
struct hcs_ident_struct {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_instance[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} hcs_ident ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident_struct))
#ifdef PROTO
int put_hcs_ident( struct hcs_ident_struct * hcs );
int get_hcs_ident( struct hcs_ident_struct * hcs );
#endif

#define NODE_IDENT_CALL 77	
#define terminate g_terminate   /* Because of name conflict */

/* end of pvm version 2.4 additions */
#endif	/* end ifdef PVM24 */


/* Added for pvm */
#ifdef PVM3
#define ANY_MSG -1              /* Any message */
#define ANY_NODE -1              /* Any node */

/* Message structures for communication between host and nodes */

/* For initialization call */
#define MAX_NUMBER_NODES 8
#define HOST_NAME_LENGTH 128
struct hcs_ident_struct {
  int msg_type;
  int node;               /* integer identifies caller's instance */
  int your_node;          /* logical node number for this node */
  int number_nodes;       /* how many nodes for this partition */
  int node_tid[MAX_NUMBER_NODES];   /* instance number for logical node */
  char host_name[HOST_NAME_LENGTH];     /* name of host */
} hcs_ident ;
#define HCS_IDENT_SIZE (sizeof(struct hcs_ident_struct))
#ifdef PROTO
int put_hcs_ident( struct hcs_ident_struct * hcs );
int get_hcs_ident( struct hcs_ident_struct * hcs );
#endif

#define NODE_IDENT_CALL 77	

/* end of pvm version 3 additions */
#endif	/* end ifdef PVM3 */

/* Added for MPL */
#ifdef MPL
#define MPL_NOT_COMPLETED -1          /* For mpc_status */
#define MPL_INACTIVE -2         /* For mpc_status */
#define DONTCARE -1             /* Any message or any node */
int mperrno;             /* Used for error reporting */
/* The following are from /usr/lpp/poe/include/mpproto.h */
extern int mpc_environ(int *howmany,int *whoami);
extern int mpc_stopall(int errcode);
extern int mpc_group(int gsize,int glist[],int label,int *gid);
extern int mpc_send(char *sarr,int len,int dest,int type,int *msgid);
extern int mpc_recv(char *darr,int len,int *src,int *type,int *msgid);
extern int mpc_bsend(char *sarr,int len,int dest,int type);
extern int mpc_brecv(char *darr,int len,int *src,int *type,int *nbytes);
extern int mpc_status(int msgid);
extern int mpc_wait(int *msgid,int *nbytes);
extern int mpc_sync(int gid);

extern void s_vadd(float in1[],float in2[],float out[],int *len);
extern void d_vadd(double in1[],double in2[],double out[],int *len);
extern void s_vmax(float in1[],float in2[],float out[],int *len);
extern void d_vmax(double in1[],double in2[],double out[],int *len);
#define MAX_NUMBER_NODES 64
#define HOST_NAME_LENGTH 128

/* end of MPL additions */
#endif	/* end ifdef MPL */

#define FIELD_REQUEST 100	/* used by field_pointer...() */
#define FIELD_REPLY 101		/* used by field_pointer...() */
#define GENERAL_GATHER_BASE 102	/* types from this to this+number_of_nodes
				are used by the general_gather routines */
#define GATHER_BASE 614	/* types greater than or equal to this are used
			   by the gather routines */
/* pid */
#define NODE_PID 0
#define HOST_PID 0
#define ALL_NODES -1	/* works for Ncube, Intel */ /* Don't use with pvm */

/* definitions of restore and save lattice commands */
#define CONTINUE 10
#define FRESH    11
#define RELOAD_ASCII  12
#define RELOAD_BINARY  13
#define RELOAD_CHECKPOINT  14
#define FORGET 20
#define SAVE_ASCII 21
#define SAVE_BINARY 22
#define SAVE_CHECKPOINT 23

/* Directions, and a macro to give the opposite direction */
/*  These must go from 0 to 7 because they will be used to index an
    array. */
/* Also define NDIRS = number of directions */
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define OPP_DIR(dir)	(7-(dir))	/* Opposite direction */
#define NDIRS 8				/* number of directions */

#define MAX_GATHERS 128	/* Maximum number of gather tables */

/* arguments to the make_gather() routine */
#define FORWARDS 1
#define BACKWARDS (-1)	/* BACKWARDS = -FORWARDS */
#define OWN_INVERSE 0
#define WANT_INVERSE 1
#define NO_INVERSE 2
#define ALLOW_EVEN_ODD 0
#define NO_EVEN_ODD 1
#define SAME_PARITY 0
#define SWITCH_PARITY 1
#define SCRAMBLE_PARITY 2

/* "comlink" is the basic structure used in gathering neighboring sites.
   Each node will maintain one such structure for each direction for each
   (other) node that contains sites that are neighbors of the sites on
   this node.  For example, if the XUP neighbors of sites on this node
   are found on two other nodes, then this node will maintain a linked
   list of two comlink structures for gathering from the XUP direction.
*/
struct comlink {
	/* pointer to next in list, NULL if this is last */
    struct comlink *nextcomlink;
	/* number of the node to which we connect */
    int othernode;
	/* number of even sites on this node that have neighbors on
	   other node connected by this "comlink", and same for odd
	   sites on this node. */
    int n_even_connected, n_odd_connected;
	/* Address of list of indices of even sites (on this node)
	   whose neighbors are found through this comlink, same for odd.
	   The odd list follows the even list, so to get all sites you
	   start at esitelist and take n_even_connected+n_odd_connected
	   addresses. */
	/* When the comlink is describing sites to be sent, the "odd"
	   list lists sites whose neighbors are even.  This convention
	   is natural for the nearest neighbor gathers.  For gathers
	   which don't allow even and odd site gathers, the even list
	   is used for list of sites to be received and the odd
	   list for sites to be sent.  Different comlink structures
	   may point to the same list.  For example, the receive list
	   for one gather may be a send list for the opposite gather. */
    int *esitelist, *ositelist;
};
typedef struct comlink comlink;


/* Structure to keep track of outstanding sends and receives */
typedef struct {
	/* node sending or receiving message */
    int msg_node;
	/* size of message in bytes */
    int msg_size;
	/* address of buffer malloc'd for message */
    char *msg_buf;
	/* message id returned by system call */
    int msg_id;
#if defined(PVM) || defined(PVM24) || defined(PVM3) || defined(MPL)
    int msg_OK;
        /* flag to track the asynchronous arrival of messages */
#endif
#ifdef MPL
    int mpl_msgid;
        /* MPL assigned message id for checking status with mpc_status */
#endif
} msg_tag;

/* Structure for requesting a field from another node */
typedef struct {
    int field;	/* offset of field in site */
    int size;	/* size of field */
    int index;	/* index of field on other node */
} msg_request;


/* Communications routines */
char * machine_type();
int mynode();
int numnodes();
void g_sync();
void g_floatsum();
void g_doublesum();
void g_floatmax();
void g_doublemax();
void g_complexsum();
void g_dcomplexsum();
void broadcast_float();
void broadcast_double();
void broadcast_complex();
void broadcast_dcomplex();
void send_integer();
void receive_integer();
double dclock();

msg_tag *start_gather();
void wait_gather();
void cleanup_gather();
msg_tag *start_general_gather();
void wait_general_gather();
void cleanup_general_gather();
char *field_pointer_at_coordinates();
char *field_pointer_at_direction();
void cleanup_field_pointer();

/* Each node maintains a list of headers to lists of comlinks */
/**EXTERN comlink * neighborlist[NDIRS];**/
/* addresses of neighboring sites, NULL if off-node */
/**EXTERN site ** neighbor[NDIRS];**/

#ifdef NONSENSE
#ifdef PE_CODE
#include <cm/cmmd.h>
/* On CM5, redefine system calls */
#define main CMPE_control
#define printf host_printf
#define scanf host_scanf
#define fprintf host_fprintf
#define fscanf host_fscanf
#define fflush host_fflush
#define fopen host_fopen
#define fclose host_fclose
#define open host_open
#define close host_close
#define creat host_creat
#define read host_read
#define write host_write
#endif
#endif

#ifdef PVM
#ifndef HOST_CODE
#define terminate g_terminate   /* Because of name conflict */
#define scanf host_scanf
#define printf host_printf
#define fflush host_fflush
#endif
#endif
