/******** layout_gen.c *********/
/* MIMD code version 3 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* generic layout - dumb */


/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
   num_sites(node) returns the number of sites on a node
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include <globaldefs.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su2.h>
#include LATDEF


/* FIRST, FOR TESTING, A COMPLETELY STUPID DISTRIBUTION */
void setup_layout(){
register int i,x,y,z,t;
#ifdef ACCORDION
    if(mynode()==0)printf("LAYOUT: Can't do ACCORDION option\n");
    exit(1);
#endif
#ifdef GRAYCODE
    if(mynode()==0)printf("LAYOUT: Can't do GRAYCODE option\n");
    exit(1);
#endif
#ifdef EVENFIRST
    if(mynode()==0)printf("LAYOUT: Can't do EVENFIRST option\n");
    exit(1);
#endif
    sites_on_node = even_sites_on_node = odd_sites_on_node = 0;
    for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++){
	if( node_number(x,y,z,t) == mynode() ){
	    sites_on_node++;
	    if( (x+y+z+t)%2 == 0) even_sites_on_node++;
	    else		  odd_sites_on_node++;
	}
    }
}
int node_number(x,y,z,t) int x,y,z,t; {
int i;
    i = x+nx*(y+ny*(z+nz*t));
    i = i%numnodes();
    return(i);
    /**return( i%numnodes() );**/
}
int node_index(x,y,z,t) int x,y,z,t; {
int i;
    i = x+nx*(y+ny*(z+nz*t));
    return( i/numnodes() );
}
int num_sites(node) int node; {
int i;
    i = nx*ny*nz*nt;
    if( node< i%numnodes() )return( i/numnodes()+1 );
    else return( i/numnodes() );
}


/* A COMPLETELY RANDOM DISTRIBUTION, FOR STRESS TESTING */
/* Only works for volume a power of 2 */
/**
void setup_layout(){
int i,x,y,z,t;
#ifdef ACCORDION
    if(mynode()==0)printf("LAYOUT: Can't do ACCORDION option\n");
    exit(1);
#endif
#ifdef GRAYCODE
    if(mynode()==0)printf("LAYOUT: Can't do GRAYCODE option\n");
    exit(1);
#endif
#ifdef EVENFIRST
    if(mynode()==0)printf("LAYOUT: Can't do EVENFIRST option\n");
    exit(1);
#endif
    i=volume; while(i%2==0)i>>=1;
    if(i!=1){printf("NEED VOLUME POWER OF TWO\n"); exit(0);}

    sites_on_node = even_sites_on_node = odd_sites_on_node = 0;
    for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++){
	if( node_number(x,y,z,t) == mynode() ){
	    sites_on_node++;
	    if( (x+y+z+t)%2 == 0) even_sites_on_node++;
	    else		  odd_sites_on_node++;
	}
    }
}
int node_number(x,y,z,t) int x,y,z,t; {
register int i,j,k;
    for(i=nx*ny*nz*nt-1,j=0;i>0;i>>=1,j++);
    for(i=numnodes()-1,k=0;i>0;i>>=1,k++);
    i = x+nx*(y+ny*(z+nz*t));
    i = (29*i+23)%(nx*ny*nz*nt);
    return( i>>(j-k) );
}
int node_index(x,y,z,t) int x,y,z,t; {
register int i,j,k;
    for(i=nx*ny*nz*nt-1,j=0;i>0;i>>=1,j++);
    for(i=numnodes()-1,k=0;i>0;i>>=1,k++);
    i = x+nx*(y+ny*(z+nz*t));
    i = (29*i+23)%(nx*ny*nz*nt);
    return( i & ((nx*ny*nz*nt-1)>>k) );
}
int num_sites(node) int node; {
register int i;
    i = nx*ny*nz*nt;
    if( node< i%numnodes() )return( i/numnodes()+1 );
    else return( i/numnodes() );
}
**/



/** A TIME SLICED DISTRIBUTION - ONE TIME SLICE TO A NODE */
/**
void setup_layout(){
#ifdef ACCORDION
    if(mynode()==0)printf("LAYOUT: Can't do ACCORDION option\n");
    exit(1);
#endif
#ifdef GRAYCODE
    if(mynode()==0)printf("LAYOUT: Can't do GRAYCODE option\n");
    exit(1);
#endif
#ifdef EVENFIRST
    if(mynode()==0)printf("LAYOUT: Can't do EVENFIRST option\n");
    exit(1);
#endif
    if(nt>numnodes()){printf("Not enough nodes\n"); terminate(1);}
    if(mynode() < nt)sites_on_node = nx*ny*nz;
    else sites_on_node =  0;
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}
int node_number(x,y,z,t) int x,y,z,t; {
    return( t );
}
int node_index(x,y,z,t) int x,y,z,t; {
    return( x+nx*(y+ny*z) );
}
int num_sites(node) int node; {
    if( node< nt )return( nx*ny*nz );
    else return( 0 );
}
**/
/** A "ZT SLICED" DISTRIBUTION: XY PLANES ALL ON ONE NODE */
/* number of slices is nz*nt, number of sites in slice is nx*ny */
/**
void setup_layout(){
}
int node_number(x,y,z,t) int x,y,z,t; {
register int i,j,k;
}
int node_index(x,y,z,t) int x,y,z,t; {
}
int num_sites(node) int node; {
    if(node < (nz*nt)%num_nodes() )return( (nx*ny)*( (nz*nt)/numnodes() +1 ) );
    else                           return( (nx*ny)*( (nz*nt)/numnodes() ) );
}
**/
