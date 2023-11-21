/******** fourier.c *************/

/* MIMD version 3.000 by T.D., D.T. */

/* These routines set up and perform Fourier transforms on fields in
   the lattice.  The field consists of "size" consecutive complex
   numbers.  For example, an su3_vector is three consecutive
   complex numbers, and a wilson_vector is 12.

   The setup_fourier() routine makes all the tables needed for the
   communications.

   References:
	Fast Fourier Transform and Convolution Algorithms, H.J. Nussbaumer
	(Springer series in information sciences, 1982)  QA 403.5
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "globaldefs.h"
#include <su2.h>
#include LATDEF /* I changed this from lattice_w.h LK */
#include <comdefs.h>  /* definitions and variables for communications */

#define FORALLUPDIR(dir) for(dir=XUP; dir<=TUP; dir++)

/* Variables set by setup_fourier() and used by fourier() */
int dim[4];		/* dimensions */
int logdim[4];		/* log_2 of dimension */
int bitrev_dir;		/* index for bit reverse gather */
int *butterfly_dir[4];	/* indices for butterflies.  First index is
			   direction, second is level.  The actual space
			   will be malloc'ed in setup_fourier(). */
			/* Level 0 is the lowest order butterfly, i.e.
			   reverse the righmost (ones) bit.  Levels
			   range from 0 to n-1 when the dimension is 2^n */

setup_fourier( key ) int *key; {
    /* "key" is a four component array.  If a component is 1, the Fourier
	transform is done in that direction, if it is 0 that direction is
	left alone. */
register int dir,i,j,k;
int arg[2];
int bitrev_map(), butterfly_map();

    /* Check that relevant dimensions are power of 2, remember their logs */
    /* Allocate space for butterfly_dir */
    dim[0]=nx;dim[1]=ny;dim[2]=nz;dim[3]=nt;
    FORALLUPDIR(dir){
	if(key[dir]==1){
	    for( j=0,i=dim[dir] ; (i&0x01) == 0; i>>=1 )j++;
/**printf("Out of loop: i=%d, j=%d\n",i,j);**/
	    if( i != 1 ){	/* Not a power of 2 */
		if(this_node==0)printf("Can't Fourier transform dimension %d\n",
		    dim[dir]);
		terminate(1);
	    }
	    logdim[dir]=j;
	    butterfly_dir[dir] = (int *)malloc( j*sizeof(int) );
	}
	else{  /* ignore this dimension */
	    logdim[dir]= -1; 
	    butterfly_dir[dir]=NULL;	/* This will remind us to ignore it */
	}
/**printf("dir = %d, dim = %d, logdim = %d\n",dir,dim[dir],logdim[dir]);**/
    }

    /* Set up bit-reverse */
    bitrev_dir = make_gather( bitrev_map, key, OWN_INVERSE,
	NO_EVEN_ODD, SCRAMBLE_PARITY );

    /* Set up butterflies */
    FORALLUPDIR(dir)if(key[dir]==1){
	for(i=0;i<logdim[dir];i++){
	    arg[0] = dir;
	    arg[1] = i;
	    butterfly_dir[dir][i] = make_gather( butterfly_map, arg,
		OWN_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY );
	}
    }
}

/* Function that defines the bit reverse mapping */
bitrev_map(x,y,z,t,key,fb,xp,yp,zp,tp)
int x,y,z,t,*key,fb,*xp,*yp,*zp,*tp;
{
   if(key[XUP]==1) *xp=bitrev_one_int(x,logdim[XUP]); else *xp=x;
   if(key[YUP]==1) *yp=bitrev_one_int(y,logdim[YUP]); else *yp=y;
   if(key[ZUP]==1) *zp=bitrev_one_int(z,logdim[ZUP]); else *zp=z;
   if(key[TUP]==1) *tp=bitrev_one_int(t,logdim[TUP]); else *tp=t;
}

/* Bit reverse a single integer, which ranges from 0 to 2^n-1 ( n bits ) */
int bitrev_one_int(i,n) int i,n; {
register int j,k;
int itemp; /*TEMPORARY*/
itemp=i; /*TEMPORARY*/
    for(k=0,j=0;j<n;j++){ /* j counts the bits, k is result we build up */
	/* This is a dumb way to do this, but it is only the setup code */
	k |= (i&0x01);	/* pull off lowest order bit of i */
	i >>= 1;	/* throw away lowest order bit of i */
	k <<= 1;	/* open up lowest order bit of j for next bit */
    }
    k >>= 1;	/* we overran by one bit */
/**printf("The %d bit reverse of %x is %x\n",n,itemp,k);**/
    return(k);
}
 
/* Function that defines the butterfly mappings */
butterfly_map(x,y,z,t,arg,fb,xp,yp,zp,tp)
int x,y,z,t,*arg,fb,*xp,*yp,*zp,*tp;
{
int i,mask,dir,level;
    /* arg contains the direction of the butterfly and the level */
    dir=arg[0];  level=arg[1];	/* just so I can remember them */
    mask = 0x01 << level;	/* one in the bit we switch */

    *xp=x; *yp=y; *zp=z; *tp=t;
    if(arg[0]==XUP) *xp ^= mask;
    if(arg[0]==YUP) *yp ^= mask;
    if(arg[0]==ZUP) *zp ^= mask;
    if(arg[0]==TUP) *tp ^= mask;
/**if(x==0 && y==0 && z==0)
printf("The level %d dir %d butterfly of %d %d %d %d is %d %d %d %d\n",
level,dir,x,y,z,t,*xp,*yp,*zp,*tp);**/
}




/* The actual Fourier transform routine.
   The algorithm is the "decimation in frequency" scheme depicted in
   Nussbaumer figure 4.2.
*/
void fourier(src,space,size,isign) 
field_offset src,space;	/* src is field to be transformed */
			/* space is working space, same size as src */
int size;		/* Size of field in bytes.  The field must
			   consist of size/sizeof(complex) consecutive
			   complex numbers.  For example, an su3_vector
			   is 3 complex numbers. */
int isign;		/* 1 for x -> k, -1 for k -> x */
{
/* local variables */
register complex *space_pt,*src_pt;
register complex cc1,cc2;
register int mask,level,i,j,n,power,dir;
register site *s;
register float theta_0;
complex *phase;	/* array of phase factors */
msg_tag *tag;
int ncomp;	/* number of complex numbers in field */

    ncomp = size/sizeof(complex);
    /* Danielson-Lanczos section */
    /* loop over all directions, and if we are supposed to transform in
	that direction, do it */
    FORALLUPDIR(dir)if(logdim[dir] != -1){
        /* The fundamental angle, others are multiples of this */
        theta_0 = -isign*2*PI/dim[dir];
	/* Make an array of phase factors */
	phase = (complex *)malloc( (dim[dir]/2)*sizeof(complex) );
	for(i=0;i<dim[dir]/2;i++)phase[i]=ce_itheta( i*theta_0 );

	for(level=logdim[dir]-1,mask=dim[dir]>>1; level>=0; level--,mask>>=1 ){
	    /* "mask" picks out the bit that is flipped to find the
		coordinate you are combining with */

	    /* Get the site at other end of butterfly */
	    tag = start_gather( src, size,
		butterfly_dir[dir][level], EVENANDODD, gen_pt[0]);
	    wait_gather(tag);
	    FORALLSITES(i,s){
		memcpy( F_PT(s,space), gen_pt[0][i], size );
	    }
	    cleanup_gather(tag);

	    FORALLSITES(i,s){
		/* Find coordinate - treat s->x,y,z,t as array */
		n = ((short *)&(s->x))[dir];
		src_pt = (complex *)F_PT(s,src);	/* pointer to source */
		space_pt = (complex *)F_PT(s,space);	/* pointer to partner */

		/* If this is the "top" site in the butterfly, just
		   add in the partner.  If it is the "bottom" site,
		   subtract site from its partner, then multiply by
		   the phase factor. */
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Local = ( %.04f , %.04f ) partner = ( %.04f , %.04f )\n",
src_pt[0].real, src_pt[0].imag, space_pt[0].real, space_pt[0].imag );
}**/
		if( n&mask ){	/* Bottom site - the bit is one */
		    if(level==0){
		        /* Special case level 0 - all phases are 1 */
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], src_pt[j] );
			}
		    }
		    else {	/* General level */
		    	power = (n&(mask-1))<<(logdim[dir]-level-1);
/**if(s->x==0 && s->y==0 && s->z==0)
printf("Level %d site %d power is %d\n", level,n,power);**/
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], cc2 );
			    /* cc2  <-  partner - local */
		            CMUL( phase[power], cc2, src_pt[j] );
			}
		    } /* end general level */
		}
		else {		/* Top site */
		    for(j=0;j<ncomp;j++){	/* loop over complex numbers */
		        CSUM( src_pt[j], space_pt[j] );
		    }
		}
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Result = ( %.04f , %.04f )\n",
src_pt[0].real, src_pt[0].imag );
}**/
	    }	/* end loop over sites */
	}  /* end loop over level */
	free(phase);
    } /* for loop on direction */

    /* Bit reverse */
    tag = start_gather( src, size, bitrev_dir,
	EVENANDODD, gen_pt[0]);
    wait_gather(tag);
    FORALLSITES(i,s){
	memcpy( F_PT(s,space), gen_pt[0][i], size );
    }
    FORALLSITES(i,s){
	memcpy( F_PT(s,src), F_PT(s,space), size );
    }
    cleanup_gather(tag);
}



#ifdef NONESENSE

/********************** TEST ROUTINE ***********************************/
test_fourier(int kz, int kt){
int i,j,k;
complex cc1,cc2;
site *s;
float phase;
int key[4];


    /* Set up Fourier transform */
    key[XUP] = 0;
    key[YUP] = 0;
    key[ZUP] = 1;
    key[TUP] = 1;
    setup_fourier(key);

    /* Set source equal to e^{ikx} */
    FORALLSITES(i,s){
        /* k = 0,0,1,1 */
	s->chi.d[0].c[0] 
           = ce_itheta( 2*PI*((float)(kt*s->t)/nt+(float)(kz*s->z)/nz) );
/**cc1 = ce_itheta(  7*2*PI*(s->t)/nt );
CSUM(s->chi.d[0].c[0],cc1);**/

	/**if(s->t==1)s->chi.d[0].c[0].real = 1.0;
	else s->chi.d[0].c[0].real = 0.0;
	s->chi.d[0].c[0].imag = 0.0;**/

/**if(s->x==0 && s->y==0 && s->z==0 )
printf("At site %d %d %d %d started with ( %.04f , %.04f )\n",
  s->x,s->y,s->z,s->t,s->chi.d[0].c[0].real,s->chi.d[0].c[0].imag);**/
    }

    /* Dump the the initial field */
    FORALLSITES(i,s){
if(s->x==0 && s->y==0 )
	printf("At site %d %d %d %d got ( %.04f , %.04f )\n",
	    s->x,s->y,s->z,s->t,s->chi.d[0].c[0].real,s->chi.d[0].c[0].imag);
    }
    printf("\n");

    /* Transform it.  Use psi as workspace */
    fourier( F_OFFSET(chi), F_OFFSET(psi), sizeof(complex), FORWARDS);

    /* Dump the result */
    FORALLSITES(i,s){
if(s->x==0 && s->y==0 )
	printf("At site %d %d %d %d got ( %.04f , %.04f )\n",
	    s->x,s->y,s->z,s->t,s->chi.d[0].c[0].real,s->chi.d[0].c[0].imag);
    }
}

#endif



