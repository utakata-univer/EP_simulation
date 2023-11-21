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

/* random number routines */
void initialize_prn();
float myrand();
typedef struct {
   /* We assume long is at least 32 bits */
    unsigned long r0,r1,r2,r3,r4,r5,r6;
    unsigned long multiplier,addend,ic_state;
   float scale;
} double_prn;


/*
/* ROUTINES */

float gaussian_rand_no(void *prn_pt);
/*       prn_pt is pointer to random # state, passed to myrand
*       file "gaussrand.c" */

