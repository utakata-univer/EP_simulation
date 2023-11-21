/*********************** io_lat.c *************************/
/* MIMD version 3 */

/* routines for lattice input/output. */
/* This works for both Intel and Ncube, but other machines may need
   special treatment */

#ifdef CMMD_2
#include <cm/cmmd.h>
#include <cm/cmmd-io.h> 
#endif	/* like this, for CM5 */
#ifdef CMMD_3
#include <cm/cmmd.h>
#endif	/* like this, for CM5 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "globaldefs.h"
#include <su2.h>
#include LATDEF
#include <comdefs.h>
#define EPS 1e-6
#include <errno.h>
extern int errno;

/* read and write an ascii lattice */
/* format:
    version_number (int)
    nx ny nz nt (int)
    coupling_1 coupling_2 (float)
    for(t=...)for(z=...)for(y=...)for(x=...){
	xlink,ylink,zlink,tlink
    }
        for each link:
            for(i=...)for(j=...){link[i][j].real, link[i][j].imag}
*/
void restore_ascii(filenam,c1,c2) char *filenam; float c1,c2; {
  FILE *fopen(),*fp;
  int currentnode,newnode;
  int version_number,i,x,y,z,t,dir;
  float x1,x2;
  su2_matrix lbuf[4];

  if(this_node==0){
    fp = fopen(filenam,"r");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filenam,errno);
      terminate(1);
    }
#ifdef CMMD_2
    CMMD_fset_io_mode(fp, CMMD_independent);
#endif
#ifdef CMMD_3
    CMMD_fset_io_mode(fp, CMMD_independent);
#endif
    if( (fscanf(fp,"%d",&version_number))!=1 ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    if(version_number != VERSION_NUMBER){
      printf("Incorrect version number in lattice header\n");terminate(1);
    }
    if( (fscanf(fp,"%d%d%d%d",&x,&y,&z,&t))!=4 ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    if( x!=nx || y!=ny || z!=nz || t!=nt ){
      printf("Incorrect lattice size %d,%d,%d,%d\n",x,y,z,t);terminate(1);
    }
    if( (fscanf(fp,"%e%e",&x1,&x2))!=2 ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    if( fabs(x1-c1) > EPS || fabs(x2-c2) > EPS ){ printf(
	  "Warning: lattice couplings, %e, %e, not equal to program's\n",
	  (double)x1,(double)x2);
    }
  }
  currentnode=0;
  g_sync();

  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    newnode=node_number(x,y,z,t);
    if(newnode != currentnode){
      g_sync();
      currentnode=newnode;
    }

/* Node 0 reads, and sends site to correct node */
    if(this_node==0){
      for(dir=XUP;dir<=TUP;dir++){
	if( fscanf(fp,"%e%e%e%e\n",
		   &(lbuf[dir].e[0]),
		   &(lbuf[dir].e[1]),
		   &(lbuf[dir].e[2]),
		   &(lbuf[dir].e[3]) )!= 4) {
	      printf("Read error in restore_ascii\n"); terminate(1);
	 }
      }
      if(currentnode==0){	/* just copy links */
	i = node_index(x,y,z,t);
        for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
      }
      else {		/* send to correct node */
	send_field(lbuf,4*sizeof(su2_matrix),currentnode);
      }
    }
  /* The node which contains this site reads message */
    else {	/* for all nodes other than node 0 */
      if(this_node==currentnode){
	get_field(lbuf,4*sizeof(su2_matrix));
	i = node_index(x,y,z,t);
	for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
      }
    }
  }

  g_sync();
  if(this_node==0){
    printf("Restored ascii lattice from file  %s\n",filenam);
    fclose(fp);
    fflush(stdout);
  }
}

void save_ascii(filenam,c1,c2) char *filenam; float c1,c2; {
FILE *fopen(),*fp;
int currentnode,newnode;
int i,j,x,y,z,t,dir;
su2_matrix lbuf[4];
    /* node 0 does all the writing */
    if(this_node==0){
        fp = fopen(filenam,"w");
        if(fp==NULL){
	    printf("Can't open file %s, error %d\n",filenam,errno);terminate(1);
        }
#ifdef CMMD_2
CMMD_fset_io_mode(fp, CMMD_independent);
#endif
#ifdef CMMD_3
CMMD_fset_io_mode(fp, CMMD_independent);
#endif
        if( (fprintf(fp,"%d\n",VERSION_NUMBER))==EOF ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
        if( (fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==EOF ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
        if( (fprintf(fp,"%.7e\t%.7e\n",(double)c1,(double)c2))==EOF ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
    }
    g_sync();
    currentnode=0;

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
	newnode=node_number(x,y,z,t);
	if(newnode != currentnode){	/* switch to another node */
	    g_sync();
	    currentnode=newnode;
	}

	if(this_node==0){
	    if(currentnode==0){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++) lbuf[dir]=lattice[i].link[dir];
	    }
	    else{
		get_field(lbuf,4*sizeof(su2_matrix));
	    }
	    for(dir=XUP;dir<=TUP;dir++){
	      if( fprintf(fp,"%e\t%e\t%e\t%e\n",
			     lbuf[dir].e[0],
			     lbuf[dir].e[1],
			     lbuf[dir].e[2],
			     lbuf[dir].e[3] ) == EOF) {
		    printf("Write error in save_ascii\n"); terminate(1);
	      }
	    }
	}
	else {	/* for nodes other than 0 */
	    if(this_node==currentnode){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
		send_field(lbuf,4*sizeof(su2_matrix),0);
	    }
	}
    }
    g_sync();
    if(this_node==0){
	fflush(fp);
	printf("Saved ascii lattice in file  %s\n",filenam);
	fclose(fp);
	fflush(stdout);
    }
}


/* read and write a binary lattice */
#include <fcntl.h>
void restore_binary(filenam,c1,c2) char *filenam; float c1,c2; {
int fd;
int currentnode,newnode;
int version_number,i,j,x,y,z,t,dir,dims[4];
float x1,x2;
su2_matrix lbuf[4];
    if(this_node==0){
        fd = open(filenam,O_RDONLY,0);
        if(fd < 0){
	    printf("Can't open file %s, error %d\n",filenam,errno);
	    terminate(1);
        }
#ifdef CMMD_2
CMMD_set_io_mode(fd, CMMD_independent);
#endif
#ifdef CMMD_3
CMMD_set_io_mode(fd, CMMD_independent);
#endif
        if( (read(fd,&version_number,sizeof(int)))!=sizeof(int) ){
	    printf("Error in reading lattice header\n"); terminate(1);
        }
        if(version_number != VERSION_NUMBER){
	    printf("Incorrect version number in lattice header\n");terminate(1);
        }
        if( (read(fd,dims,4*sizeof(int)))!=4*sizeof(int) ){
	    printf("Error in reading lattice header\n"); terminate(1);
        }
	x = dims[XUP]; y=dims[YUP]; z=dims[ZUP]; t=dims[TUP];
        if( x!=nx || y!=ny || z!=nz || t!=nt ){
	    printf("Incorrect lattice size %d,%d,%d,%d\n",x,y,z,t);terminate(1);
        }
        if( (read(fd,&x1,sizeof(float)))!=sizeof(float) ){
	    printf("Error in reading lattice header\n"); terminate(1);
        }
        if( (read(fd,&x2,sizeof(float)))!=sizeof(float) ){
	    printf("Error in reading lattice header\n"); terminate(1);
        }
        if( fabs(x1-c1) > EPS || fabs(x2-c2) > EPS ){ printf(
            "Warning: lattice couplings, %e, %e, not equal to program's\n",
                (double)x1,(double)x2);
        }
    }
    currentnode=0;
    g_sync();

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
	newnode=node_number(x,y,z,t);
	if(newnode != currentnode){
	    g_sync();
	    currentnode=newnode;
	}

	/* Node 0 reads, and sends site to correct node */
	if(this_node==0){
	    if( (read(fd,lbuf,4*sizeof(su2_matrix))) != 4*sizeof(su2_matrix)){
		printf("Read error in restore_binary\n"); terminate(1);
	    }
	    if(currentnode==0){	/* just copy links */
		i = node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
	    }
	    else {		/* send to correct node */
		send_field(lbuf,4*sizeof(su2_matrix),currentnode);
	    }
	}

	/* The node which contains this site reads message */
	else {	/* for all nodes other than node 0 */
	    if(this_node==currentnode){
		get_field(lbuf,4*sizeof(su2_matrix));
		i = node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
	    }
	}
    }

    g_sync();
    if(this_node==0){
	printf("Restored binary lattice from file  %s\n",filenam);
        close(fd);
	fflush(stdout);
    }
}

void save_binary(filenam,c1,c2) char *filenam; float c1,c2; {
int fd;
int currentnode,newnode;
int i,j,x,y,z,t,dir,dims[4];
su2_matrix lbuf[4];
float x1,x2;
    /* node 0 does all the writing */
    if(this_node==0){
        fd = creat(filenam,0644);
        if(fd < 0){
	    printf("Can't open file %s, error %d\n",filenam,errno);terminate(1);
        }
#ifdef CMMD_2
CMMD_set_io_mode(fd, CMMD_independent);
#endif
#ifdef CMMD_3
CMMD_set_io_mode(fd, CMMD_independent);
#endif
	i=VERSION_NUMBER;
        if( (write(fd,&i,sizeof(int))) != sizeof(int) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	dims[XUP]=nx; dims[YUP]=ny; dims[ZUP]=nz; dims[TUP]=nt;
        if( (write(fd,dims,4*sizeof(int))) != 4*sizeof(int) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	x1 = c1; x2 = c2; /* can't take address of argument */
        if( (write(fd,&x1,sizeof(float))) != sizeof(float) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
        if( (write(fd,&x2,sizeof(float))) != sizeof(float) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
    }
    g_sync();
    currentnode=0;

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
	newnode=node_number(x,y,z,t);
	if(newnode != currentnode){	/* switch to another node */
	    g_sync();
	    currentnode=newnode;
	}

	if(this_node==0){
	    if(currentnode==0){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
	    }
	    else{
		get_field(lbuf,4*sizeof(su2_matrix));
	    }
	    if( (write(fd,lbuf,4*sizeof(su2_matrix))) != 4*sizeof(su2_matrix) ){
		printf("Write error in save_binary\n"); terminate(1);
	    }
	}
	else {	/* for nodes other than 0 */
	    if(this_node==currentnode){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
		send_field(lbuf,4*sizeof(su2_matrix),0);
	    }
	}
    }
    g_sync();
    if(this_node==0){
	printf("Saved binary lattice in file  %s\n",filenam);
	close(fd);
	fflush(stdout);
    }
}

/* read and write a checkpoint lattice */
/* This is not implemented in general, but is used for machines
   with special, usually parallel, file systems */
#include <fcntl.h>
void restore_checkpoint(filenam,c1,c2) char *filenam; float c1,c2; {
int fd;
int currentnode,newnode;
int version_number,i,j,x,y,z,t,dims[4];
float x1,x2;
su2_matrix lbuf[4];

    if(this_node==0){
	printf("restore_checkpoint is not implemented\n"); terminate(1);
    }
}

void save_checkpoint(filenam,c1,c2) char *filenam; float c1,c2; {
int fd;
int currentnode,newnode;
int i,j,x,y,z,t,dir,dims[4];
su2_matrix lbuf[4];
char *writebuf;
float x1,x2;
double tsize;

    if(this_node==0){
	printf("save_checkpoint is not implemented\n"); terminate(1);
    }
}

