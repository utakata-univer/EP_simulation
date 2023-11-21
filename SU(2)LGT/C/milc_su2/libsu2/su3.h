/******************************  su3.h **********************************
*									*
*  Defines and subroutine declarations for SU3 simulation		*
*  MIMD version 3 							*
*									*
*/
typedef struct { complex e[3][3]; } su3_matrix;
typedef struct { complex c[3]; } su3_vector;
typedef struct
  { complex m01,m02,m12; float m00im,m11im,m22im; float space; } anti_hermitmat;
typedef struct { su3_vector d[4]; } wilson_vector;
typedef struct { su3_vector h[2]; } half_wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;

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
* ROUTINES FOR SU(3) MATRIX OPERATIONS
*
* void mult_su3_nn( a,b,c )
*	su3_matrix *a,*b,*c;
*	matrix multiply, no adjoints
*	files "m_mat_nn.c", "m_mat_nn.m4"
* void mult_su3_na( a,b,c )
*	su3_matrix *a,*b,*c;
*	matrix multiply, second matrix is adjoint
*	files "m_mat_na.c", "m_mat_na.m4"
* void mult_su3_an( a,b,c )
*	su3_matrix *a,*b,*c;
*	matrix multiply, first matrix is adjoint
*	files "m_mat_an.c", "m_mat_an.m4"
* float realtrace_su3(a,b)
*	su3_matrix *a,*b;  (Re(Tr( A_adjoint*B)) )
*	file "realtr.c"
* complex trace_su3(a)
*	su3_matrix *a; 
*	file "trace_su3.c"
* complex det_su3(a)
*	su3_matrix *a;
*	file "det_su3.c"
* void add_su3_matrix(a,b,c)
*	su3_matrix *a,*b,*c;
*	file "addmat.c"
* void sub_su3_matrix(a,b,c)
*	su3_matrix *a,*b,*c;
*	file "submat.c"
* void scalar_mult_add_su3_matrix(a,b,s,c)
*	su3_matrix *a,*b,*c; float s;
*	file "s_m_a_mat.c"
* void scalar_mult_sub_su3_matrix(a,b,s,c)
*	su3_matrix *a,*b,*c; float s;
*	file "s_m_s_mat.c"
* void c_scalar_mult_add_su3mat(m1,m2,phase,m3)
*	su3_matrix *m1,*m2,*m3; complex *phase;
*	file "cs_m_a_mat.c"
* void su3_adjoint(a,b)
*	su3_matrix *a,*b;
*	file "su3_adjoint.c"
* void make_anti_hermitian(m3,ah3)
*	su3_matrix *m3; anti_hermitmat *ah3;
*	file "make_ahmat.c"
* void random_anti_hermitian(mat_antihermit,prn_pt)
*	anti_hermitmat *mat_antihermit;
*	void *prn_pt;   (passed through to myrand())
*	file "rand_ahmat.c"
* void uncompress_anti_hermitian(mat_antihermit,mat_su3)
*	anti_hermitmat *mat_antihermit; su3_matrix *mat_su3;
*	file "uncmp_ahmat.c"
* void compress_anti_hermitian(mat_su3,mat_antihermit)
*	anti_hermitmat *mat_antihermit; su3_matrix *mat_su3;
*	file "cmp_ahmat.c"
* void su3mat_copy(a,b)
*	su3_matrix *a,*b;
*	file "su3mat_copy.c"
*
*
* ROUTINES FOR su3_vector OPERATIONS ( 3 COMPONENT COMPLEX )
*
* void c_scalar_mult_su3vec(v1,phase,v2)
*	su3_vector *v1,*v2; complex *phase;
*	file "cs_m_vec.c"
* void c_scalar_mult_add_su3vec(v1,phase,v2)
*	su3_vector *v1,*v2; complex *phase;
*	file "cs_m_a_vec.c"
* void c_scalar_mult_sub_su3vec(v1,phase,v2)
*	su3_vector *v1,*v2; complex *phase;
*	file "cs_m_s_vec.c"
* void su3_projector(a,b,c)
*	su3_vector *a,*b; su3_matrix *c;
*	( outer product of A and B)
*	file "su3_proj.c"
* void su3vec_copy(a,b)
*	su3_vector *a,*b;
*	file "su3vec_copy.c"
* 
* void mult_su3_mat_vec( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_matvec.c", "m_matvec.m4"
* void mult_su3_mat_vec_sum( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_matvec_s.c", "m_matvec_s.m4"
* void mult_su3_mat_vec_sum_4dir( a,b0,b1,b2,b3,c )
*	su3_matrix *a; su3_vector *b0,*b1,*b2,*b3,*c;
*	file "m_mv_s_4dir.c", "m_mv_s_4dir.m4"
*	file "m_mv_s_4di2.m4" is alternate version with pipelined loads.
*	Multiply four su3_vectors by elements of an array of su3_matrices,
*	sum results.
*	C <- A[0]*B0 + A[1]*B1 + A[2]*B2 + A[3]*B3
* void mult_su3_mat_vec_nsum( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_matvec_ns.c"
* void mult_adj_su3_mat_vec( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_amatvec.c", "m_amatvec.m4"
* void mult_adj_su3_mat_vec_4dir( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_amv_4dir.c", "m_amv_4dir.m4"
*	file "m_amv_4di2.m4" is alternate version with pipelined loads.
*	Multiply an su3_vector by adjoints of elements of an array 
*	of su3_matrices, results in an array of su3_vectors.
*	C[i] <- A_adjoint[i]*B, i = 0,1,2,3
* void mult_adj_su3_mat_vec_sum( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_amatvec_s.c"
* void mult_adj_su3_mat_vec_nsum( a,b,c )
*	su3_matrix *a; su3_vector *b,*c;
*	file "m_amatvec_ns.c"
* void add_su3_vector(a,b,c)
*	su3_vector *a,*b,*c;
*	file "addvec.c", "addvec.m4"
* void sub_su3_vector(a,b,c)
*	su3_vector *a,*b,*c;
*	file "subvec.c", "subvec.m4"
* void sub_four_su3_vecs(a,b1,b2,b3,b4)
*	su3_vector *a,*b1,*b2,*b3,*b4;
*	file "sub4vecs.c", "sub4vecs.m4"
* void scalar_mult_su3_vector(a,s,c)
*	su3_vector *a,*c; float s;
*	file "s_m_vec.c"
* void scalar_mult_add_su3_vector(a,b,s,c)
*	su3_vector *a,*b,*c; float s;
*	file "s_m_a_vec.c", "s_m_a_vec.m4"
* void scalar_mult_sum_su3_vector(a,b,s)
*	su3_vector *a,*b; float s;
*	file "s_m_s_vec.c", "s_m_s_vec.m4"
* void scalar_mult_sub_su3_vector(a,b,s,c)
*	su3_vector *a,*b,*c; float s;
*	file "s_m_s_vec.c"
* complex su3_dot(a,b)
*	su3_vector *a,*b;
*	file "su3_dot.c"
* float su3_rdot(a,b)
*	su3_vector *a,*b;
*	file "su3_rdot.c", "su3_rdot.m4"
* float magsq_su3vec(a)
*	su3_vector *a;
*	file "msq_su3vec.c", "msq_su3vec.m4"
*
* ROUTINES FOR WILSON VECTORS
*
* void mult_mat_wilson_vec();
*	file m_mat_wvec.c
*	   mult_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
*		wilson_vector *dest)
*	   dest <- mat*src
* void mult_su3_mat_hwvec();
*	file m_mat_hwvec.c
*	   mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
*		half_wilson_vector *dest)
*	   dest <- mat*src
* void mult_adj_mat_wilson_vec();
*	file m_amat_wvec.c
*	   mult_adj_mat_wilson_vec( su3_matrix *mat, wilson_vector *src,
*		wilson_vector *dest)
*	   dest <- mat_adjoint*src
* void mult_adj_su3_mat_hwvec();
*	file m_amat_hwvec.c
*	   mult_adj_su3_mat_hwvec( su3_matrix *mat,
*		half_wilson_vector *src, half_wilson_vector *dest)
*	   dest <- mat_adjoint*src
*
* void add_wilson_vector();
*	file add_wvec.c
*	   add_wilson_vector( wilson_vector *src1,*src2,*dest)
*	   dest <- src1+src2
* void sub_wilson_vector();
*	file sub_wvec.c
*	   sub_wilson_vector( wilson_vector *src1,*src2,*dest)
*	   dest <- src1-src2
* void scalar_mult_wvec();
*	file s_m_wvec.c
*	   scalar_mult_wvec( wilson_vector *src, float s, wilson_vector *dest)
*	   dest <- s*src
* void scalar_mult_hwvec();
*	file s_m_hwvec.c
*	   scalar_mult_hwvec( half_wilson_vector *src, float s,
*		half_wilson_vector *dest)
*	   dest <- s*src
* float magsq_wvec();
*	file msq_wvec.c
*	   s = magsq_wvec( wilson_vector *src)
*	   s <- squared magnitude of src
* complex wvec_dot();
*	file wvec_dot.c
*	   c = wvec_dot( wilson_vector *src1, *src2 )
*	   c <- dot product of src1 and src2
* float wvec_rdot(a,b)
*	wilson_vector *a,*b;
*	file "wvec_rdot.c", "wvec_rdot.m4"
*	   r = wvec_dot( wilson_vector *src1, *src2 )
*	   r <- real part of dot product of src1 and src2
* void scalar_mult_add_wvec();
*	file s_m_a_wvec.c
*	   scalar_mult_add_wvec( wilson_vector *src1,*src2, float s,
*		wilson_vector *dest)
*	   dest <- src1 + s*src2
* void scalar_mult_addtm_wvec();
*	file s_m_atm_wvec.c
*	   scalar_mult_addtm_wvec( wilson_vector *src1,*src2, float s,
*		wilson_vector *dest)
*	   dest <- -src1 + s*src2   ("atm"="add to minus")
* void c_scalar_mult_add_wvec(v1,v2,phase,v3)
*       wilson_vector *v1,*v2,*v3; complex *phase;
*       file "cs_m_a_wvec.c"
* void su3_projector_w(a,b,c)
*	wilson_vector *a,*b; su3_matrix *c;
*	sum over spins of outer product of A.d[s] and B.d[s]  - a three
*	  by three complex matrix
*	file "su3_proj_w.c"
* void copy_wvec();
*	file copy_wvec.c
*	   copy_wvec( wilson_vector *src, wilson_vector *dest)
*	   dest <- src
* void clear_wvec();
*	file clear_wvec.c
*	  clear_wvec( wilson_vector *dest)
*	   dest <- 0.0
* void wp_shrink();
*	file wp_shrink.c , wp_shrink.m4
*	   wp_shrink( wilson_vector *src, half_wilson_vector *dest,
*		int dir, int sign)
*	   if(dir = [XYZT]UP) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
* void wp_shrink_4dir();
*	file wp_shrink4.c wp_shrink4.m4
*	  void wp_shrink_4dir(a,b1,b2,b3,b4,sign)
*		wilson_vector *a;
*		half_wilson_vector *b1,*b2,*b3,*b4;
*		int sign;
*	  Shrink A in X,Y,Z,T directions respectively, results in B1-B4
* void wp_grow();
*	file wp_grow.c , wp_grow.m4
*	   wp_grow( half_wilson_vector *src, wilson_vector *dest,
*		int dir, int sign)
*	   if(dir = [XYZT]UP) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
*  	    Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to
*		multiplication by 1+-gamma_dir, or 1-+gamma_dir if sign=MINUS
* void wp_grow_add();
*	file wp_grow_a.c , wp_grow_a.m4
*	   wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
*		int dir, int sign)
*	   wp_grow, and add result to previous contents of dest.
* void grow_add_four_wvecs();
*	file grow4wvecs.c grow4wvecs.m4
*	  void grow_add_four_wvecs(a,b1,b2,b3,b4,sign,sum)
*		wilson_vector *a;
*		half_wilson_vector *b1,*b2,*b3,*b4;
*		int sign,sum;
*         If sum==0
*	  Grow b1-b4 in X,Y,Z,T directions respectively, sum of results in A
*         If sum==1
*	  Grow b1-b4 in X,Y,Z,T directions respectively, add to current A
* void mult_by_gamma();
*	file mb_gamma.c
*	   mult_by_gamma( wilson_vector *src,*dest, int dir )
*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
* void mult_by_gamma_left();
*	file mb_gamma_l.c
*	   mult_by_gamma_left( wilson_matrix *src,*dest, int dir )
*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
*	   acts on first index of matrix
* void dump_wilson_vec();
*	file dump_wvec.c
*	   dump_wilson_vec( wilson_vector *src);
*	   print out a wilson vector
*
* MISCELLANEOUS ROUTINES
*
* float gaussian_rand_no(prn_pt)
*	void *prn_pt;  ( passed to myrand())
*	file "gaussrand.c"
*
* void dumpmat(m)
*	su3_matrix *m;
*	file "dumpmat.c"
* void dumpvec(v)
*	su3_vector *v;
*	file "dumpvec.c"
*/

void mult_su3_nn ();
void mult_su3_na ();
void mult_su3_an ();
float realtrace_su3();
complex trace_su3();
complex det_su3();
void add_su3_matrix();
void sub_su3_matrix();
void su3_adjoint();
void make_anti_hermitian();
void random_anti_hermitian();
void uncompress_anti_hermitian();
void compress_anti_hermitian();
void su3mat_copy();
void dumpmat();
void c_scalar_mult_add_su3mat();

void c_scalar_mult_su3vec    ();
void c_scalar_mult_add_su3vec();
void c_scalar_mult_sub_su3vec();
void su3_projector();
void su3vec_copy();
void mult_su3_mat_vec();
void mult_su3_mat_vec_sum();
void mult_su3_mat_vec_sum_4dir();
void mult_su3_mat_vec_nsum();
void mult_adj_su3_mat_vec();
void mult_adj_su3_mat_vec_4dir();
void mult_adj_su3_mat_vec_sum();
void mult_adj_su3_mat_vec_nsum();
void add_su3_vector();
void sub_su3_vector();
void sub_four_su3_vecs();
void scalar_mult_su3_vector();
#ifdef PROTO
void scalar_mult_add_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_sum_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar);
void scalar_mult_sub_su3_vector( su3_vector *src1, su3_vector *src2,
	float scalar, su3_vector *dest);
void scalar_mult_add_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
void scalar_mult_sub_su3_matrix( su3_matrix *src1, su3_matrix *src2,
	float scalar, su3_matrix *dest);
void scalar_mult_add_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest);
void scalar_mult_addtm_wvec( wilson_vector *src1, wilson_vector *src2,
	float scalar, wilson_vector *dest);
void c_scalar_mult_add_su3vec(su3_vector *v1, complex *phase, su3_vector
*v2);
void c_scalar_mult_add_wvec(wilson_vector *src1, wilson_vector *src2, complex *phase, wilson_vector *dest);
#else
void scalar_mult_add_su3_vector();
void scalar_mult_sum_su3_vector();
void scalar_mult_sub_su3_vector();
void scalar_mult_add_su3_matrix();
void scalar_mult_sub_su3_matrix();
void scalar_mult_add_wvec();
void scalar_mult_addtm_wvec();
void c_scalar_mult_add_su3vec();
void c_scalar_mult_add_wvec();
#endif
complex su3_dot();
float su3_rdot();
float magsq_su3vec();
void dumpvec();

void mult_mat_wilson_vec();
void mult_su3_mat_hwvec();
void mult_adj_mat_wilson_vec();
void mult_adj_su3_mat_hwvec();
void dump_wilson_vec();
void add_wilson_vector();
void scalar_mult_wvec();
void scalar_mult_hwvec();
float magsq_wvec();
complex wvec_dot();
float wvec_rdot();
void clear_wvec();
void wp_shrink();
void wp_shrink_4dir();
void wp_grow();
void wp_grow_add();
void grow_add_four_wvecs();
void mult_by_gamma();
void mult_by_gamma_left();
void su3_projector_w();
void copy_wvec();

float gaussian_rand_no();

