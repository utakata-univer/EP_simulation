/***********************  globaldefs.h **********************************
*									*
*  Defines and subroutines common to several field theories		*
*  MIMD version 4 							*
*									*
*/

#define GAMMAFIVE -1    /* some integer which is not a direction */
#define PLUS 1          /* flags for selecting M or M_adjoint */
#define MINUS -1
/* Macros to multiply complex numbers by +-1 and +-i */
#define TIMESPLUSONE(a,b) { (b).real =  (a).real; (b).imag = (a).imag; }
#define TIMESMINUSONE(a,b) { (b).real =  -(a).real; (b).imag = -(a).imag; }
#define TIMESPLUSI(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
#define TIMESMINUSI(a,b) { (b).real =  (a).imag; (b).imag = -(a).real; }

#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/
#define CMINUSEQ(a,b) (a).real -= (b).real; (a).imag -= (b).imag;
/* a += b */

#define CPLUSEQ(a,b)  (a).real += (b).real; (a).imag += (b).imag;
/* a -= b */

#define CMAGSQ(z) (((z).real)*((z).real) + ((z).imag)*((z).imag))

#define PRINTSITECOORDS(s) printf("site [%d,%d,%d,%d]:\n", \
                                     s->x,s->y,s->z,s->t);
#define PRINTSITECOORDSN(s) printf("site [%d,%d,%d,%d]: ", \
                                     s->x,s->y,s->z,s->t);



/* random number routines */
void initialize_prn();
double myrand();
typedef struct {
   /* We assume long is at least 32 bits */
    unsigned long r0,r1,r2,r3,r4,r5,r6;
    unsigned long multiplier,addend,ic_state;
   double scale;
} double_prn;


/*
/* ROUTINES */

double gaussian_rand_no(void *prn_pt);
/*       prn_pt is pointer to random # state, passed to myrand
*       file "gaussrand.c" */

