/******************************  su2.h **********************************
*									*
*  Defines and subroutine declarations for SU2 simulation		*
*  MIMD version 4,   J. Hetrick and D. Toussaint			*
*									*
* We assume ANSI c compilers, so we can use prototypes			*
* You must also include globaldefs.h					*
*									*
*/

#ifndef _SU2_H
#define _SU2_H

typedef struct { float e[4]; } su2_matrix;
     /* su2_matrix = e[0]$1 + i*e[1]*\tau_1 ... */
     /* really only "SU(2)" if a[0]^2+a[1]^2+a[2]^2+a[3]^2=1 */
     /* For (anti)Hermitian matrices, have a[0]=0 */
typedef su2_matrix su2_anti_hermitmat;
typedef struct { complex c[2]; } su2_vector;
typedef struct { su2_vector d[4]; } su2_wilson_vector;
typedef struct { su2_vector h[2]; } su2_half_wilson_vector;
typedef struct { su2_wilson_vector c[2]; } su2_color_wilson_vector;
typedef struct { su2_color_wilson_vector d[4]; } su2_wilson_matrix;


/*-----------------------------------------------------------------
* ROUTINES FOR SU(2) MATRIX OPERATIONS
* -----------------------------------------------------------------*/

void mult_su2_nn(su2_matrix *a, su2_matrix *b, su2_matrix *c );
/*	matrix multiply, no adjoints
*	files "m_mat_nn.c" */

void mult_su2_na(su2_matrix *a, su2_matrix *b, su2_matrix *c);
/*	matrix multiply, second matrix is adjoint
*	files "m_mat_na.c" */

void mult_su2_an(su2_matrix *a, su2_matrix *b, su2_matrix *c);
/*	matrix multiply, first matrix is adjoint
*	files "m_mat_an.c" */

void mult_su2_aa(su2_matrix *a, su2_matrix *b, su2_matrix *c);
/*	matrix multiply, both adjoint
*	files "m_mat_aa.c" */

void mult_su2_choose(su2_matrix *a, int flag_a, su2_matrix *b, int flag_b,
                     su2_matrix *c);
/*	matrix multiply, if flag_[ab] is -1, corresponding matrix is adjoint
*	(set flag to +1 for non-adjoint)
*	files "m_mat_choose.c" */

float realtrace_su2(su2_matrix *a, su2_matrix *b);
/*	file "realtr.c" */

float trace_su2(su2_matrix *a);
/*	su2_matrix *a; 
*	file "trace_su2.c" */

float det_su2(su2_matrix *a);
/*	file "det_su2.c" */

void add_su2_matrix(su2_matrix *a, su2_matrix *b, su2_matrix *c);
/*	file "addmat.c" */

void sub_su2_matrix(su2_matrix *a, su2_matrix *b, su2_matrix *c);
/*	file "submat.c" */

void scalar_mult_add_su2_matrix(su2_matrix *a, su2_matrix *b, float s,
                                su2_matrix *c);
/*	file "s_m_a_mat.c" */

void scalar_mult_sub_su2_matrix(su2_matrix *a, su2_matrix *b, float s,
                                su2_matrix *c);
/*	file "s_m_s_mat.c" */

void c_scalar_mult_add_su2mat(su2_matrix *a, su2_matrix *b, complex phase,
                              su2_matrix *c);
/*	file "cs_m_a_mat.c" */

void su2_adjoint(su2_matrix *a, su2_matrix *b);
/*	file "su2_adjoint.c" */

void make_anti_hermitian(su2_matrix *m2, su2_anti_hermitmat *ah2);
/*	file "make_ahmat.c"     --- takes traceless part */

void compress_anti_hermitian(su2_matrix *mat_su2, su2_anti_hermitmat *ah2);
/*	file "cmp_ahmat.c"  - anithermitian part, including trace
*	--- nothing to do */

void random_anti_hermitian(su2_anti_hermitmat *mat_antihermit, void *prn_pt);
/*	void *prn_pt;   (passed through to myrand())
*	file "rand_ahmat.c" */

void su2mat_copy(su2_matrix *a, su2_matrix *b);
/*	file "su2mat_copy.c" */

void reunit_su2(su2_matrix *a);
/*      reunitarize su2_matrix so that det a = 1 */
/*      file "reunit_su2.c"   */

su2_matrix make_su2_matrix(float a0, float a1, float a2, float a3);
/*      returns the su2_matrix made from a0..a3  */
/*      file "make_su2_mat.c"  */


/* -----------------------------------------------------------------*/
/* ROUTINES FOR su2_vector OPERATIONS ( 2 COMPONENT COMPLEX ) */
/* -----------------------------------------------------------------*/

void c_scalar_mult_su2vec(su2_vector *v1, complex *phase, su2_vector *v2);
/*	file "cs_m_vec.c" */

void c_scalar_mult_add_su2vec(su2_vector *v1, complex *phase, su2_vector *v2);
/*	file "cs_m_a_vec.c" */

void c_scalar_mult_sub_su2vec(su2_vector *v1, complex *phase, su2_vector *v2);
/*	file "cs_m_s_vec.c" */


/* We have no data type which holds the outer product of A ^ B (\in SL2(C)?)
* Thus this funtion doesn't work yet. Fortunately we don't need it.
* void su2_projector(su2_vector *a, su2_vector *b, su2_matrix *c);
*	outer product of A and B
*	file "su2_proj.c" 
* Thus we use the following:  */
void su2_antiherm_projector(su2_vector *a, su2_vector *b, su2_matrix *c);
/*      ( antihermitian part of outer product of A and B)
*       file "su2_anti_proj.c" */

void su2vec_copy(su2_vector *a, su2_vector *b);
/*	file "su2vec_copy.c" */

void mult_su2_mat_vec(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_matvec.c" */

void mult_su2_mat_vec_sum(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_matvec_s.c" */

void mult_su2_mat_vec_sum_4dir(su2_matrix *a, su2_vector *b0, su2_vector *b1,
                               su2_vector *b2, su2_vector *b3, su2_vector *c);
/*	file "m_mv_s_4dir.c", 
*	Multiply four su2_vectors by elements of an array of su2_matrices,
*	sum results.
*	C <- A[0]*B0 + A[1]*B1 + A[2]*B2 + A[3]*B3 */

void mult_su2_mat_vec_nsum(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_matvec_ns.c" */

void mult_adj_su2_mat_vec(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_amatvec.c" */

void mult_adj_su2_mat_vec_4dir(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_amv_4dir.c" 
*	Multiply an su2_vector by adjoints of elements of an array 
*	of su2_matrices, results in an array of su2_vectors.
*	C[i] <- A_adjoint[i]*B, i = 0,1,2,3 */

void mult_adj_su2_mat_vec_sum(su2_matrix *a, su2_vector *b, su2_vector *c);
/*	file "m_amatvec_s.c" */

void mult_adj_su2_mat_vec_nsum(su2_matrix *a, su2_vector *b,su2_vector *c);
/*	file "m_amatvec_ns.c" */

void add_su2_vector(su2_vector *a, su2_vector *b, su2_vector *c);
/*	file "addvec.c" */

void sub_su2_vector(su2_vector *a, su2_vector *b, su2_vector *c);
/*	file "subvec.c" */

void sub_four_su2_vecs(su2_vector *a, su2_vector *b1, su2_vector *b2,
                       su2_vector *b3, su2_vector *b4);
/*	file "sub4vecs.c" */

void scalar_mult_su2_vector(su2_vector *a, float s, su2_vector *c);
/*	file "s_m_vec.c" */

void scalar_mult_add_su2_vector(su2_vector *a, su2_vector *b, float s,
				su2_vector *c);
/*	file "s_m_a_vec.c" */

void scalar_mult_sub_su2_vector(su2_vector *a, su2_vector *b, float s, 
				su2_vector *c);
/*	file "s_m_s_vec.c" */

void scalar_mult_sum_su2_vector(su2_vector *a, su2_vector *b, float s);
/*	file "s_m_sum_vec.c" */

complex su2_dot(su2_vector *a, su2_vector *b);
/*	file "su2_dot.c" */

float su2_rdot(su2_vector *a, su2_vector *b);
/*	file "su2_rdot.c" */

float magsq_su2vec(su2_vector *a);
/*	file "msq_su2vec.c" */




/*---------------------------implimented above here------------*/




/* ROUTINES FOR WILSON VECTORS */

void mult_mat_su2_wilson_vec(su2_matrix *mat, su2_wilson_vector *src,
			     su2_wilson_vector *dest);
/*	   dest <- mat*src
*	file m_mat_wvec.c */


void mult_mat_su2_hwvec(su2_matrix *mat, su2_half_wilson_vector *src,
			su2_half_wilson_vector *dest);
/*	   dest <- mat*src
*	file m_mat_hwvec.c */

void mult_adj_mat_su2_wilson_vec( su2_matrix *mat, su2_wilson_vector *src,
				 su2_wilson_vector *dest);
/*	   dest <- mat_adjoint*src
*	file m_amat_wvec.c */

void mult_adj_mat_su2_hwvec(su2_matrix *mat, su2_half_wilson_vector *src, 
			    su2_half_wilson_vector *dest);
/*	   dest <- mat_adjoint*src
*	file m_amat_hwvec.c */

void add_su2_wilson_vector(su2_wilson_vector *src1,
			   su2_wilson_vector *src2,
			   su2_wilson_vector *dest);
/*	   dest <- src1+src2
*	file add_wvec.c */

void sub_su2_wilson_vector(su2_wilson_vector *src1,
			   su2_wilson_vector *src2, 
			   su2_wilson_vector *dest);
/*	   dest <- src1-src2
*	file sub_wvec.c */

void scalar_mult_su2_wvec(su2_wilson_vector *src, float s,
			  su2_wilson_vector *dest);
/*	   dest <- s*src
*	file s_m_wvec.c */

void scalar_mult_su2_hwvec( su2_half_wilson_vector *src, float s,
			   su2_half_wilson_vector *dest);
/*	   dest <- s*src
*	file s_m_hwvec.c */

float magsq_su2_wvec(su2_wilson_vector *src);
/*	   s <- squared magnitude of src
*	file msq_wvec.c */

complex su2_wvec_dot(su2_wilson_vector *src1, su2_wilson_vector *src2);
/*	   c <- dot product of src1 and src2
*	file wvec_dot.c */

float su2_wvec_rdot(su2_wilson_vector *src1, su2_wilson_vector *src2);
/*	   r <- real part of dot product of src1 and src2
*	file "wvec_rdot.c" */

void scalar_mult_add_su2_wvec(su2_wilson_vector *src1,
			      su2_wilson_vector *src2, float s,
			      su2_wilson_vector *dest);
/*	   dest <- src1 + s*src2
*	file s_m_a_wvec.c
 */


void scalar_mult_addtm_su2_wvec(su2_wilson_vector *src1,
				su2_wilson_vector *src2, float s,
				su2_wilson_vector *dest);
/*	   dest <- -src1 + s*src2   ("atm"="add to minus")
*	file s_m_atm_wvec.c */

void c_scalar_mult_add_su2_wvec(su2_wilson_vector *v1, 
				su2_wilson_vector *v2,
				complex *phase, su2_wilson_vector *v3);
/*       file "cs_m_a_wvec.c" */


void su2_antiherm_projector_w(su2_wilson_vector *a, su2_wilson_vector *b,
		     su2_matrix *c);
/*	sum over spins of outer product of A.d[s] and B.d[s]  - a two
*	  by two complex matrix
* This version projects onto the traceless antihermitian part.
*	file "su2_antproj_w.c" */

void copy_su2_wvec(su2_wilson_vector *src, su2_wilson_vector *dest);
/*	   dest <- src
*	file copy_wvec.c */

void clear_su2_wvec(su2_wilson_vector *dest);
/*	   dest <- 0.0
*	file clear_wvec.c */

void su2_wp_shrink(su2_wilson_vector *src, su2_half_wilson_vector *dest,
		   int dir, int sign);
/*	   if(dir = [XYZT]UP) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src along eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
*	file wp_shrink.c */


void su2_wp_shrink_4dir(su2_wilson_vector *a,
			su2_half_wilson_vector *b1, 
			su2_half_wilson_vector *b2, 
			su2_half_wilson_vector *b3, 
			su2_half_wilson_vector *b4, int sign);
/*	  Shrink A in X,Y,Z,T directions respectively, results in B1-B4
*	file wp_shrink4.c  */


void su2_wp_grow(su2_half_wilson_vector *src, su2_wilson_vector *dest,
		 int dir, int sign);
/*	   if(dir = [XYZT]UP) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue +1
*	   if(dir = [XYZT]DOWN) dest <- components of src times eigenvectors
*		of gamma_dir with eigenvalue -1
*	   if(sign==MINUS)switch roles of +1 and -1
*  	    Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to
*		multiplication by 1+-gamma_dir, or 1-+gamma_dir if sign=MINUS
*	file wp_grow.c */ 

void su2_wp_grow_add(su2_half_wilson_vector *src,
		     su2_wilson_vector *dest, int dir, int sign);
/*	   wp_grow, and add result to previous contents of dest.
*	file wp_grow_a.c */

void grow_add_four_su2_wvecs(su2_wilson_vector *a,
			     su2_half_wilson_vector *b1,
			     su2_half_wilson_vector *b2,
			     su2_half_wilson_vector *b3,
			     su2_half_wilson_vector *b4, int sign, int sum);
/*         If sum==0
*	  Grow b1-b4 in X,Y,Z,T directions respectively, sum of results in A
*         If sum==1
*	  Grow b1-b4 in X,Y,Z,T directions respectively, add to current A
*	file grow4wvecs.c */


void mult_su2_by_gamma(su2_wilson_vector *src, su2_wilson_vector *dest,
		       int dir );
/*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
*	file mb_gamma.c */


void mult_su2_by_gamma_left(su2_wilson_matrix *src, su2_wilson_matrix *dest,
			    int dir);
/*	   dest <- gamma[dir] * src,  dir=[XYZT]UP,GAMMAFIVE
*	   acts on first index of matrix
*	file mb_gamma_l.c */

void mult_su2_by_gamma_right(su2_wilson_matrix *src, su2_wilson_matrix *dest,
			    int dir);
/*	   dest <- src * gamma[dir],  dir=[XYZT]UP,GAMMAFIVE
*	   acts on second index of matrix
*	file mb_gamma_r.c */


void dump_su2_wilson_vec(su2_wilson_vector *src);
/*	   print out a wilson vector
*	file dump_wvec.c */


/* MISCELLANEOUS ROUTINES */

void dump_su2_mat(su2_matrix *m);
/*	file "dumpmat.c" */

void dump_su2_vec(su2_vector *v);
/*	file "dumpvec.c" */

#endif /* _SU2_H */
