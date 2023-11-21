
!  Program qcdf90, module field_algebra, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

!  This module defines overloaded operator that perform arithmetic
!  operations between fields and other variables.
!  The simble usead are:
!  g gauge_field
!  u full_gauge_field
!  f fermi_field
!  c complex_field
!  r real_field
!  ge generator_field 
!  m 2X2 matrix
!  complex complex variable of kind REAL8
!  real real variable of kind REAL8
!
!  integer loops indices are: xyzt for the lattice sites,
!  i,j for color indicies, s for spin index.
!


MODULE field_algebra

   USE precisions
   USE constants
   USE global_module
   IMPLICIT NONE



   INTERFACE OPERATOR(+)
      MODULE PROCEDURE g_plus_g, g_plus_m, m_plus_g, f_plus_f, &
                       c_plus_c, c_plus_r, r_plus_c, c_plus_complex, &
                       complex_plus_c, c_plus_real, real_plus_c, &
                       r_plus_r, r_plus_real, real_plus_r, &
                       ge_plus_ge, m_plus_m
   END INTERFACE

   INTERFACE OPERATOR(-)
      MODULE PROCEDURE g_minus_g, g_minus_m, m_minus_g, f_minus_f, &
                       c_minus_c, c_minus_r, r_minus_c, c_minus_complex, &
                       complex_minus_c, c_minus_real, real_minus_c, &
                       r_minus_r, r_minus_real, real_minus_r, &
                       ge_minus_ge, m_minus_m
   END INTERFACE

   INTERFACE OPERATOR(*)
      MODULE PROCEDURE  g_times_g, g_times_f, f_times_g, g_times_c, &
                        c_times_g, g_times_r, r_times_g, g_times_m, &
                        m_times_g, g_times_complex, complex_times_g, &
                        g_times_real, real_times_g, f_times_c, &
                        c_times_f, f_times_r, r_times_f, &
                        f_times_m, m_times_f, f_times_complex, &
                        complex_times_f, f_times_real, real_times_f, &
                        c_times_c, c_times_r, &
                        r_times_c, c_times_complex, complex_times_c, &
                        c_times_real, real_times_c, r_times_r, &
                        r_times_real, real_times_r, ge_times_r, &
                        r_times_ge, ge_times_real, real_times_ge, &
                        m_times_m, m_times_complex, complex_times_m, &
                        m_times_real, real_times_m, &
                        f_contract_f, ge_contract_ge
   END INTERFACE

   INTERFACE OPERATOR(/)
      MODULE PROCEDURE  g_divided_c, g_divided_r, g_divided_complex, &
                        g_divided_real, f_divided_c, f_divided_r, &
                        f_divided_complex, f_divided_real, &
                        c_divided_c, c_divided_r, r_divided_c, &
                        c_divided_complex, complex_divided_c, &
                        c_divided_real, real_divided_c,r_divided_r, &
                        r_divided_real, real_divided_r, &
                        ge_divided_r, ge_divided_real, &
                        m_divided_complex, m_divided_real, &
                        g_divided_g, m_divided_g, g_divided_m, &
                        f_divided_g, f_divided_m
   END INTERFACE

   INTERFACE OPERATOR(//)
      MODULE PROCEDURE  g_predivided_g, g_predivided_m, &
                        m_predivided_g, g_predivided_f, &
                        m_predivided_f, f_diadic_f 
   END INTERFACE

   INTERFACE OPERATOR(.Dot.)
      MODULE PROCEDURE g_dot_g
   END INTERFACE

   INTERFACE OPERATOR(.Gamma.)
!
! GAMMA1: 0 0 0 1   GAMMA2: 0 0 0-I   GAMMA3: 0 0 1 0   GAMMA4: 1 0 0 0
!         0 0 1 0           0 0 I 0           0 0 0-1           0 1 0 0
!         0 1 0 0           0-I 0 0           1 0 0 0           0 0-1 0
!         1 0 0 0           I 0 0 0           0-1 0 0           0 0 0-1
!
!
! GAMMA5: 0 0-I 0    GAMMA5=GAMMA1*GAMMA2*GAMMA3*GAMMA4
!         0 0 0-I    GAMMA(5k)=GAMMA5*GAMMA(k)
!         I 0 0 0    GAMMA(k5)=GAMMA(k)*GAMMA5
!         0 I 0 0    GAMMA(kj)=(I/2)(GAMMA(k)*GAMMA(j)-GAMMA(j)*GAMMA(k))
!
      MODULE PROCEDURE gamma_f, f_gamma
   END INTERFACE

  INTERFACE OPERATOR(.Lambda.)
!
! LAMBDA1: 0 1    LAMBDA2: 0-I    LAMBDA3: 1 0    
!          1 0             I 0             0-1    

      MODULE PROCEDURE  lambda_g, g_lambda
   END INTERFACE

   INTERFACE OPERATOR(.I.)
      MODULE PROCEDURE I_times_g, I_times_f, I_times_c, I_times_r
   END INTERFACE

  INTERFACE OPERATOR(.Minus.)
      MODULE PROCEDURE  minus_g, minus_f, minus_c, &
                        minus_r, minus_ge
   END INTERFACE

  INTERFACE OPERATOR(.Conjg.)
      MODULE PROCEDURE  conjg_g, conjg_f, conjg_c, conjg_m
   END INTERFACE

  INTERFACE OPERATOR(.Adj.)
      MODULE PROCEDURE  adj_g, adj_m
   END INTERFACE

  INTERFACE OPERATOR(.Ctr.)
      MODULE PROCEDURE  ctr_g, ctr_m
   END INTERFACE

  INTERFACE OPERATOR(.Tr.)
      MODULE PROCEDURE  tr_g, tr_m
   END INTERFACE

   INTERFACE OPERATOR(.Sqrt.)
      MODULE PROCEDURE sqrt_r
   END INTERFACE

   INTERFACE OPERATOR(.Exp.)
      MODULE PROCEDURE exp_r
   END INTERFACE

   INTERFACE OPERATOR(.Log.)
      MODULE PROCEDURE log_r
   END INTERFACE

   INTERFACE OPERATOR(.Cos.)
      MODULE PROCEDURE cos_r
   END INTERFACE

   INTERFACE OPERATOR(.Sin.)
      MODULE PROCEDURE sin_r
   END INTERFACE


CONTAINS

!+
   FUNCTION g_plus_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_plus_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_plus_g%parity=a%parity
      ELSE
        g_plus_g%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        g_plus_g%dir=a%dir
      ELSE
        g_plus_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_plus_g%fc(i,j,xyzt)=a%fc(i,j,xyzt)+b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_plus_g

   FUNCTION g_plus_m(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(gauge_field) g_plus_m
      INTEGER xyzt,i,j
      g_plus_m%parity=a%parity
      g_plus_m%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_plus_m%fc(i,j,xyzt)=a%fc(i,j,xyzt)+b%mc(i,j)
      END DO
      END DO
      END DO
   END FUNCTION g_plus_m

   FUNCTION m_plus_g(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) m_plus_g
      INTEGER xyzt,i,j
      m_plus_g%parity=b%parity
      m_plus_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        m_plus_g%fc(i,j,xyzt)=a%mc(i,j)+b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION m_plus_g

   FUNCTION f_plus_f(a,b)
      TYPE(fermi_field), INTENT(IN) :: a,b 
      TYPE(fermi_field) f_plus_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_plus_f%parity=a%parity
      ELSE
        f_plus_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_plus_f%fc(i,xyzt,s)=a%fc(i,xyzt,s)+b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION f_plus_f

   FUNCTION c_plus_c(a,b)
      TYPE(complex_field), INTENT(IN) :: a,b
      TYPE(complex_field) c_plus_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_plus_c%parity=a%parity
      ELSE
        c_plus_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_plus_c%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
      END DO
   END FUNCTION c_plus_c

   FUNCTION c_plus_r(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(complex_field) c_plus_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_plus_r%parity=a%parity
      ELSE
        c_plus_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_plus_r%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
      END DO
   END FUNCTION c_plus_r

   FUNCTION r_plus_c(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) r_plus_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_plus_c%parity=a%parity
      ELSE
        r_plus_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_plus_c%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
      END DO
   END FUNCTION r_plus_c

   FUNCTION c_plus_complex(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_plus_complex
      INTEGER xyzt
      c_plus_complex%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_plus_complex%fc(xyzt)=a%fc(xyzt)+b         
      END DO
   END FUNCTION c_plus_complex

   FUNCTION complex_plus_c(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) complex_plus_c
      INTEGER xyzt
      complex_plus_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        complex_plus_c%fc(xyzt)=a+b%fc(xyzt)         
      END DO
   END FUNCTION complex_plus_c

   FUNCTION c_plus_real(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_plus_real
      INTEGER xyzt
      c_plus_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_plus_real%fc(xyzt)=a%fc(xyzt)+b
      END DO
   END FUNCTION c_plus_real

   FUNCTION real_plus_c(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) real_plus_c
      INTEGER xyzt
      real_plus_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_plus_c%fc(xyzt)=a+b%fc(xyzt)         
      END DO
   END FUNCTION real_plus_c

   FUNCTION r_plus_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      TYPE(real_field) r_plus_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_plus_r%parity=a%parity
      ELSE
        r_plus_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_plus_r%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
      END DO
   END FUNCTION r_plus_r

   FUNCTION r_plus_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(real_field) r_plus_real
      INTEGER xyzt
      r_plus_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        r_plus_real%fc(xyzt)=a%fc(xyzt)+b
      END DO
   END FUNCTION r_plus_real

   FUNCTION real_plus_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(real_field) real_plus_r
      INTEGER xyzt
      real_plus_r%parity=b%parity
      DO xyzt=0,NXYZT2-1
         real_plus_r%fc(xyzt)=a+b%fc(xyzt)
      END DO
   END FUNCTION real_plus_r

   FUNCTION ge_plus_ge(a,b)
      TYPE(generator_field), INTENT(IN) :: a,b
      TYPE(generator_field) ge_plus_ge
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        ge_plus_ge%parity=a%parity
      ELSE
        ge_plus_ge%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        ge_plus_ge%dir=a%dir
      ELSE
        ge_plus_ge%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_plus_ge%fc(i,xyzt)=a%fc(i,xyzt)+b%fc(i,xyzt)
      END DO
      END DO
   END FUNCTION ge_plus_ge  

   FUNCTION m_plus_m(a,b)
      TYPE(matrix), INTENT(IN) :: a,b
      TYPE(matrix) m_plus_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_plus_m%mc(i,j)=a%mc(i,j)+b%mc(i,j)
      END DO
      END DO
   END FUNCTION m_plus_m

!-
   FUNCTION g_minus_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_minus_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_minus_g%parity=a%parity
      ELSE
        g_minus_g%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        g_minus_g%dir=a%dir
      ELSE
        g_minus_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_minus_g%fc(i,j,xyzt)=a%fc(i,j,xyzt)-b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_minus_g

   FUNCTION g_minus_m(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(gauge_field) g_minus_m
      INTEGER xyzt,i,j
      g_minus_m%parity=a%parity
      g_minus_m%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_minus_m%fc(i,j,xyzt)=a%fc(i,j,xyzt)-b%mc(i,j)
      END DO
      END DO
      END DO
   END FUNCTION g_minus_m

   FUNCTION m_minus_g(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) m_minus_g
      INTEGER xyzt,i,j
      m_minus_g%parity=b%parity
      m_minus_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        m_minus_g%fc(i,j,xyzt)=a%mc(i,j)-b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION m_minus_g

   FUNCTION f_minus_f(a,b)
      TYPE(fermi_field), INTENT(IN) :: a,b 
      TYPE(fermi_field) f_minus_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_minus_f%parity=a%parity
      ELSE
        f_minus_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_minus_f%fc(i,xyzt,s)=a%fc(i,xyzt,s)-b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION f_minus_f   

   FUNCTION c_minus_c(a,b)
      TYPE(complex_field), INTENT(IN) :: a,b
      TYPE(complex_field) c_minus_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_minus_c%parity=a%parity
      ELSE
        c_minus_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_minus_c%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
      END DO
   END FUNCTION c_minus_c  
 
   FUNCTION c_minus_r(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(complex_field) c_minus_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_minus_r%parity=a%parity
      ELSE
        c_minus_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_minus_r%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
      END DO
   END FUNCTION c_minus_r

   FUNCTION r_minus_c(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) r_minus_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_minus_c%parity=a%parity
      ELSE
        r_minus_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_minus_c%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
      END DO
   END FUNCTION r_minus_c

   FUNCTION c_minus_complex(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_minus_complex
      INTEGER xyzt
      c_minus_complex%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_minus_complex%fc(xyzt)=a%fc(xyzt)-b         
      END DO
   END FUNCTION c_minus_complex

   FUNCTION complex_minus_c(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) complex_minus_c
      INTEGER xyzt
      complex_minus_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        complex_minus_c%fc(xyzt)=a-b%fc(xyzt)
      END DO
   END FUNCTION complex_minus_c

   FUNCTION c_minus_real(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_minus_real
      INTEGER xyzt
      c_minus_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_minus_real%fc(xyzt)=a%fc(xyzt)-b
      END DO
   END FUNCTION c_minus_real

   FUNCTION real_minus_c(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) real_minus_c
      INTEGER xyzt
      real_minus_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_minus_c%fc(xyzt)=a-b%fc(xyzt)      
      END DO
   END FUNCTION real_minus_c

   FUNCTION r_minus_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      TYPE(real_field) r_minus_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_minus_r%parity=a%parity
      ELSE
        r_minus_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_minus_r%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
      END DO
   END FUNCTION r_minus_r

   FUNCTION r_minus_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(real_field) r_minus_real
      INTEGER xyzt
      r_minus_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        r_minus_real%fc(xyzt)=a%fc(xyzt)-b
      END DO
   END FUNCTION r_minus_real

   FUNCTION real_minus_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(real_field) real_minus_r
      INTEGER xyzt
      real_minus_r%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_minus_r%fc(xyzt)=a-b%fc(xyzt)
      END DO
   END FUNCTION real_minus_r

   FUNCTION ge_minus_ge(a,b)
      TYPE(generator_field), INTENT(IN) :: a,b
      TYPE(generator_field) ge_minus_ge
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        ge_minus_ge%parity=a%parity
      ELSE
        ge_minus_ge%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        ge_minus_ge%dir=a%dir
      ELSE
        ge_minus_ge%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_minus_ge%fc(i,xyzt)=a%fc(i,xyzt)-b%fc(i,xyzt)
      END DO
      END DO
   END FUNCTION ge_minus_ge

   FUNCTION m_minus_m(a,b)
      TYPE(matrix), INTENT(IN) :: a,b
      TYPE(matrix) m_minus_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_minus_m%mc(i,j)=a%mc(i,j)-b%mc(i,j)
      END DO
      END DO
   END FUNCTION m_minus_m

!*
   FUNCTION g_times_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_times_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_times_g%parity=a%parity
      ELSE
        g_times_g%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        g_times_g%dir=a%dir
      ELSE
        g_times_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_g%fc(i,j,xyzt)= &
                       a%fc(i,1,xyzt)*b%fc(1,j,xyzt)+ &
                       a%fc(i,2,xyzt)*b%fc(2,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_times_g

   FUNCTION g_times_f(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) g_times_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        g_times_f%parity=a%parity
      ELSE
        g_times_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        g_times_f%fc(i,xyzt,s)= &
                       a%fc(i,1,xyzt)*b%fc(1,xyzt,s)+ &
                       a%fc(i,2,xyzt)*b%fc(2,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION g_times_f

   FUNCTION f_times_g(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(fermi_field) f_times_g
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_times_g%parity=a%parity
      ELSE
        f_times_g%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_g%fc(i,xyzt,s)= &
                       a%fc(1,xyzt,s)*b%fc(1,i,xyzt)+ &
                       a%fc(2,xyzt,s)*b%fc(2,i,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION f_times_g

   FUNCTION g_times_c(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(gauge_field) g_times_c
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_times_c%parity=a%parity
      ELSE
        g_times_c%parity=-1
      ENDIF
      g_times_c%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_c%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_times_c

   FUNCTION c_times_g(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) c_times_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        c_times_g%parity=a%parity
      ELSE
        c_times_g%parity=-1
      ENDIF
      c_times_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        c_times_g%fc(i,j,xyzt)=a%fc(xyzt)*b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION c_times_g

   FUNCTION g_times_r(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(gauge_field) g_times_r
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_times_r%parity=a%parity
      ELSE
        g_times_r%parity=-1
      ENDIF
      g_times_r%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_r%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_times_r

   FUNCTION r_times_g(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) r_times_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        r_times_g%parity=a%parity
      ELSE
        r_times_g%parity=-1
      ENDIF
      r_times_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        r_times_g%fc(i,j,xyzt)=a%fc(xyzt)*b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION r_times_g

   FUNCTION g_times_m(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(gauge_field) g_times_m
      INTEGER xyzt,i,j
      g_times_m%parity=a%parity
      g_times_m%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_m%fc(i,j,xyzt)=a%fc(i,1,xyzt)*b%mc(1,j)+ &
                               a%fc(i,2,xyzt)*b%mc(2,j)
      END DO
      END DO
      END DO
   END FUNCTION g_times_m

   FUNCTION m_times_g(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) m_times_g
      INTEGER xyzt,i,j
      m_times_g%parity=b%parity
      m_times_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        m_times_g%fc(i,j,xyzt)=a%mc(i,1)*b%fc(1,j,xyzt)+ &
                               a%mc(i,2)*b%fc(2,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION m_times_g

   FUNCTION g_times_complex(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(gauge_field) :: g_times_complex
      INTEGER :: xyzt,i,j
      g_times_complex%parity=a%parity
      g_times_complex%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_complex%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b
      END DO
      END DO
      END DO
   END FUNCTION g_times_complex

   FUNCTION complex_times_g(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) :: complex_times_g
      INTEGER :: xyzt,i,j
      complex_times_g%parity=b%parity
      complex_times_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        complex_times_g%fc(i,j,xyzt)=a*b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION complex_times_g

   FUNCTION g_times_real(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(gauge_field) :: g_times_real
      INTEGER :: xyzt,i,j
      g_times_real%parity=a%parity
      g_times_real%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_times_real%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b
      END DO
      END DO
      END DO
   END FUNCTION g_times_real

   FUNCTION real_times_g(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) :: real_times_g
      INTEGER :: xyzt,i,j
      real_times_g%parity=b%parity
      real_times_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        real_times_g%fc(i,j,xyzt)=a*b%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION real_times_g

   FUNCTION f_times_c(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(fermi_field) f_times_c
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_times_c%parity=a%parity
      ELSE
        f_times_c%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_c%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION f_times_c

   FUNCTION c_times_f(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) c_times_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        c_times_f%parity=a%parity
      ELSE
        c_times_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        c_times_f%fc(i,xyzt,s)=a%fc(xyzt)*b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION c_times_f

   FUNCTION f_times_r(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(fermi_field) f_times_r
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_times_r%parity=a%parity
      ELSE
        f_times_r%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_r%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION f_times_r

   FUNCTION r_times_f(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) r_times_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        r_times_f%parity=a%parity
      ELSE
        r_times_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        r_times_f%fc(i,xyzt,s)=a%fc(xyzt)*b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION r_times_f

   FUNCTION f_times_m(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(fermi_field) f_times_m
      INTEGER s,xyzt,i
      f_times_m%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_m%fc(i,xyzt,s)=a%fc(1,xyzt,s)*b%mc(1,i)+ &
                               a%fc(2,xyzt,s)*b%mc(2,i)
      END DO
      END DO
      END DO
   END FUNCTION f_times_m

   FUNCTION m_times_f(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) m_times_f
      INTEGER s,xyzt,i
      m_times_f%parity=b%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        m_times_f%fc(i,xyzt,s)=a%mc(i,1)*b%fc(1,xyzt,s)+ &
                               a%mc(i,2)*b%fc(2,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION m_times_f

   FUNCTION f_times_complex(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(fermi_field) f_times_complex
      INTEGER s,xyzt,i
      f_times_complex%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_complex%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b
      END DO
      END DO
      END DO
   END FUNCTION f_times_complex

   FUNCTION complex_times_f(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) complex_times_f
      INTEGER s,xyzt,i
      complex_times_f%parity=b%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        complex_times_f%fc(i,xyzt,s)=a*b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION complex_times_f

   FUNCTION f_times_real(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(fermi_field) f_times_real
      INTEGER s,xyzt,i
      f_times_real%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_times_real%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b
      END DO
      END DO
      END DO
   END FUNCTION f_times_real

   FUNCTION real_times_f(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) real_times_f
      INTEGER s,xyzt,i
      real_times_f%parity=b%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        real_times_f%fc(i,xyzt,s)=a*b%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION real_times_f

   FUNCTION c_times_c(a,b)
      TYPE(complex_field), INTENT(IN) :: a,b
      TYPE(complex_field) c_times_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_times_c%parity=a%parity
      ELSE
        c_times_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_times_c%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
      END DO
   END FUNCTION c_times_c

   FUNCTION c_times_r(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(complex_field) c_times_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_times_r%parity=a%parity
      ELSE
        c_times_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_times_r%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
      END DO
   END FUNCTION c_times_r

   FUNCTION r_times_c(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) r_times_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_times_c%parity=a%parity
      ELSE
        r_times_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_times_c%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
      END DO
   END FUNCTION r_times_c   

   FUNCTION c_times_complex(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_times_complex
      INTEGER xyzt
      c_times_complex%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_times_complex%fc(xyzt)=a%fc(xyzt)*b
      END DO
   END FUNCTION c_times_complex

   FUNCTION complex_times_c(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) complex_times_c
      INTEGER xyzt
      complex_times_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        complex_times_c%fc(xyzt)=a*b%fc(xyzt)
      END DO
   END FUNCTION complex_times_c

   FUNCTION c_times_real(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_times_real
      INTEGER xyzt
      c_times_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_times_real%fc(xyzt)=a%fc(xyzt)*b
      END DO
   END FUNCTION c_times_real

   FUNCTION real_times_c(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) real_times_c
      INTEGER xyzt
      real_times_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_times_c%fc(xyzt)=a*b%fc(xyzt)
      END DO
   END FUNCTION real_times_c

   FUNCTION r_times_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      TYPE(real_field) r_times_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_times_r%parity=a%parity
      ELSE
        r_times_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_times_r%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
      END DO
   END FUNCTION r_times_r

   FUNCTION r_times_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(real_field) r_times_real
      INTEGER xyzt
      r_times_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        r_times_real%fc(xyzt)=a%fc(xyzt)*b
      END DO
   END FUNCTION r_times_real

   FUNCTION real_times_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(real_field) real_times_r
      INTEGER xyzt
      real_times_r%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_times_r%fc(xyzt)=a*b%fc(xyzt)
      END DO
   END FUNCTION real_times_r

   FUNCTION ge_times_r(a,b)
      TYPE(generator_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(generator_field) ge_times_r
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        ge_times_r%parity=a%parity
      ELSE
        ge_times_r%parity=-1
      ENDIF
      ge_times_r%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_times_r%fc(i,xyzt)=a%fc(i,xyzt)*b%fc(xyzt)
      END DO
      END DO
   END FUNCTION ge_times_r

   FUNCTION r_times_ge(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(generator_field), INTENT(IN) :: b
      TYPE(generator_field) r_times_ge
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        r_times_ge%parity=a%parity
      ELSE
        r_times_ge%parity=-1
      ENDIF
      r_times_ge%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        r_times_ge%fc(i,xyzt)=a%fc(xyzt)*b%fc(i,xyzt)
      END DO
      END DO
   END FUNCTION r_times_ge

   FUNCTION ge_times_real(a,b)
      TYPE(generator_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(generator_field) :: ge_times_real
      INTEGER xyzt,i
      ge_times_real%parity=a%parity
      ge_times_real%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_times_real%fc(i,xyzt)=a%fc(i,xyzt)*b
      END DO
      END DO
   END FUNCTION ge_times_real

   FUNCTION real_times_ge(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(generator_field), INTENT(IN) :: b
      TYPE(generator_field) :: real_times_ge
      INTEGER xyzt,i
      real_times_ge%parity=b%parity
      real_times_ge%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        real_times_ge%fc(i,xyzt)=a*b%fc(i,xyzt)
      END DO
      END DO
   END FUNCTION real_times_ge

   FUNCTION m_times_m(a,b)
      TYPE(matrix), INTENT(IN) :: a,b
      TYPE(matrix) m_times_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_times_m%mc(i,j)=a%mc(i,1)*b%mc(1,j)+ &
                          a%mc(i,2)*b%mc(2,j)
      END DO
      END DO
   END FUNCTION m_times_m

   FUNCTION m_times_complex(a,b)
      TYPE(matrix), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(matrix) m_times_complex
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_times_complex%mc(i,j)=a%mc(i,j)*b
      END DO
      END DO
   END FUNCTION m_times_complex

   FUNCTION complex_times_m(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(matrix) complex_times_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        complex_times_m%mc(i,j)=a*b%mc(i,j)
      END DO
      END DO
   END FUNCTION complex_times_m

   FUNCTION m_times_real(a,b)
      TYPE(matrix), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(matrix) m_times_real
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_times_real%mc(i,j)=a%mc(i,j)*b
      END DO
      END DO
   END FUNCTION m_times_real

   FUNCTION real_times_m(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(matrix) real_times_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        real_times_m%mc(i,j)=a*b%mc(i,j)
      END DO
      END DO
   END FUNCTION real_times_m

   FUNCTION f_contract_f(a,b)
      TYPE(fermi_field), INTENT(IN) :: a,b
      TYPE(complex_field) f_contract_f
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        f_contract_f%parity=a%parity
      ELSE
        f_contract_f%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
         f_contract_f%fc(xyzt)=CONJG(a%fc(1,xyzt,1))*b%fc(1,xyzt,1) &
                           +CONJG(a%fc(2,xyzt,1))*b%fc(2,xyzt,1) &
                           +CONJG(a%fc(1,xyzt,2))*b%fc(1,xyzt,2) &
                           +CONJG(a%fc(2,xyzt,2))*b%fc(2,xyzt,2) &
                           +CONJG(a%fc(1,xyzt,3))*b%fc(1,xyzt,3) &
                           +CONJG(a%fc(2,xyzt,3))*b%fc(2,xyzt,3) &
                           +CONJG(a%fc(1,xyzt,4))*b%fc(1,xyzt,4) &
                           +CONJG(a%fc(2,xyzt,4))*b%fc(2,xyzt,4) 
      END DO
   END FUNCTION f_contract_f

   FUNCTION ge_contract_ge(a,b)
      TYPE(generator_field), INTENT(IN) :: a,b
      TYPE(real_field) ge_contract_ge
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        ge_contract_ge%parity=a%parity
      ELSE
        ge_contract_ge%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        ge_contract_ge%fc(xyzt)= &
                     a%fc(1,xyzt)*b%fc(1,xyzt)
        DO i=2,3
          ge_contract_ge%fc(xyzt)=ge_contract_ge%fc(xyzt)+ &
                     a%fc(i,xyzt)*b%fc(i,xyzt)
        END DO
      END DO
   END FUNCTION ge_contract_ge

!/
   FUNCTION g_divided_c(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(gauge_field) g_divided_c
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_divided_c%parity=a%parity
      ELSE
        g_divided_c%parity=-1
      ENDIF
      g_divided_c%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_c%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_divided_c

   FUNCTION g_divided_r(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(gauge_field) g_divided_r
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_divided_r%parity=a%parity
      ELSE
        g_divided_r%parity=-1
      ENDIF
      g_divided_r%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_r%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_divided_r

   FUNCTION g_divided_complex(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(gauge_field) :: g_divided_complex
      INTEGER :: xyzt,i,j
      g_divided_complex%parity=a%parity
      g_divided_complex%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_complex%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b
      END DO
      END DO
      END DO
   END FUNCTION g_divided_complex

   FUNCTION g_divided_real(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(gauge_field) :: g_divided_real
      INTEGER :: xyzt,i,j
      g_divided_real%parity=a%parity
      g_divided_real%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_real%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b
      END DO
      END DO
      END DO
   END FUNCTION g_divided_real

   FUNCTION f_divided_c(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_c
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_divided_c%parity=a%parity
      ELSE
        f_divided_c%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_c%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION f_divided_c

   FUNCTION f_divided_r(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_r
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_divided_r%parity=a%parity
      ELSE
        f_divided_r%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_r%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b%fc(xyzt)
      END DO
      END DO
      END DO
   END FUNCTION f_divided_r

   FUNCTION f_divided_complex(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_complex
      INTEGER s,xyzt,i
      f_divided_complex%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_complex%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b
      END DO
      END DO
      END DO
   END FUNCTION f_divided_complex

   FUNCTION f_divided_real(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_real
      INTEGER s,xyzt,i
      f_divided_real%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_real%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b
      END DO
      END DO
      END DO
   END FUNCTION f_divided_real

   FUNCTION c_divided_c(a,b)
      TYPE(complex_field), INTENT(IN) :: a,b
      TYPE(complex_field) c_divided_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_divided_c%parity=a%parity
      ELSE
        c_divided_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_divided_c%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
      END DO
   END FUNCTION c_divided_c

   FUNCTION c_divided_r(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(complex_field) c_divided_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        c_divided_r%parity=a%parity
      ELSE
        c_divided_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        c_divided_r%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
      END DO
   END FUNCTION c_divided_r

   FUNCTION r_divided_c(a,b)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) r_divided_c
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_divided_c%parity=a%parity
      ELSE
        r_divided_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_divided_c%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
      END DO
   END FUNCTION r_divided_c

   FUNCTION c_divided_complex(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_divided_complex
      INTEGER xyzt
      c_divided_complex%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_divided_complex%fc(xyzt)=a%fc(xyzt)/b
      END DO
   END FUNCTION c_divided_complex

   FUNCTION complex_divided_c(a,b)
      COMPLEX(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) complex_divided_c
      INTEGER xyzt
      complex_divided_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        complex_divided_c%fc(xyzt)=a/b%fc(xyzt)
      END DO
   END FUNCTION complex_divided_c

   FUNCTION c_divided_real(a,b)
      TYPE(complex_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(complex_field) c_divided_real
      INTEGER xyzt
      c_divided_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        c_divided_real%fc(xyzt)=a%fc(xyzt)/b
      END DO
   END FUNCTION c_divided_real

   FUNCTION real_divided_c(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) real_divided_c
      INTEGER xyzt
      real_divided_c%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_divided_c%fc(xyzt)=a/b%fc(xyzt)
      END DO
   END FUNCTION real_divided_c

   FUNCTION r_divided_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      TYPE(real_field) r_divided_r
      INTEGER xyzt
      IF(a%parity.EQ.b%parity) THEN
        r_divided_r%parity=a%parity
      ELSE
        r_divided_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        r_divided_r%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
      END DO
   END FUNCTION r_divided_r

   FUNCTION r_divided_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(real_field) r_divided_real
      INTEGER xyzt
      r_divided_real%parity=a%parity
      DO xyzt=0,NXYZT2-1
        r_divided_real%fc(xyzt)=a%fc(xyzt)/b
      END DO
   END FUNCTION r_divided_real

   FUNCTION real_divided_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(real_field) real_divided_r
      INTEGER xyzt
      real_divided_r%parity=b%parity
      DO xyzt=0,NXYZT2-1
        real_divided_r%fc(xyzt)=a/b%fc(xyzt)
      END DO
   END FUNCTION real_divided_r

   FUNCTION ge_divided_r(a,b)
      TYPE(generator_field), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(generator_field) ge_divided_r
      INTEGER xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        ge_divided_r%parity=a%parity
      ELSE
        ge_divided_r%parity=-1
      ENDIF
      ge_divided_r%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_divided_r%fc(i,xyzt)=a%fc(i,xyzt)/b%fc(xyzt)
      END DO
      END DO
   END FUNCTION ge_divided_r

   FUNCTION ge_divided_real(a,b)
      TYPE(generator_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(generator_field) :: ge_divided_real
      INTEGER :: xyzt,i
      ge_divided_real%parity=a%parity
      ge_divided_real%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        ge_divided_real%fc(i,xyzt)=a%fc(i,xyzt)/b
      END DO
      END DO
   END FUNCTION ge_divided_real

   FUNCTION m_divided_complex(a,b)
      TYPE(matrix), INTENT(IN) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      TYPE(matrix) m_divided_complex
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_divided_complex%mc(i,j)=a%mc(i,j)/b
      END DO
      END DO
   END FUNCTION m_divided_complex

   FUNCTION m_divided_real(a,b)
      TYPE(matrix), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      TYPE(matrix) m_divided_real
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        m_divided_real%mc(i,j)=a%mc(i,j)/b
      END DO
      END DO
   END FUNCTION m_divided_real

   FUNCTION g_divided_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_divided_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_divided_g%parity=a%parity
      ELSE
        g_divided_g%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        g_divided_g%dir=a%dir
      ELSE
        g_divided_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_g%fc(i,j,xyzt)= &
            a%fc(i,1,xyzt)*CONJG(b%fc(j,1,xyzt))+ &
            a%fc(i,2,xyzt)*CONJG(b%fc(j,2,xyzt))
      END DO
      END DO
      END DO
   END FUNCTION g_divided_g

   FUNCTION m_divided_g(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) m_divided_g
      INTEGER xyzt,i,j
      m_divided_g%parity=b%parity
      m_divided_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        m_divided_g%fc(i,j,xyzt)= &
                 a%mc(i,1)*CONJG(b%fc(j,1,xyzt))+ &
                 a%mc(i,2)*CONJG(b%fc(j,2,xyzt))
      END DO
      END DO
      END DO
   END FUNCTION m_divided_g

   FUNCTION g_divided_m(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(gauge_field) g_divided_m
      INTEGER xyzt,i,j
      g_divided_m%parity=a%parity
      g_divided_m%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_divided_m%fc(i,j,xyzt)= &
                 a%fc(i,1,xyzt)*CONJG(b%mc(j,1))+ &
                 a%fc(i,2,xyzt)*CONJG(b%mc(j,2))
      END DO
      END DO
      END DO
   END FUNCTION g_divided_m

   FUNCTION f_divided_g(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_g
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        f_divided_g%parity=a%parity
      ELSE
        f_divided_g%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_g%fc(i,xyzt,s)= &
            a%fc(1,xyzt,s)*CONJG(b%fc(i,1,xyzt))+ &
            a%fc(2,xyzt,s)*CONJG(b%fc(i,2,xyzt))
      END DO
      END DO
      END DO
   END FUNCTION f_divided_g

   FUNCTION f_divided_m(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(fermi_field) f_divided_m
      INTEGER s,xyzt,i
      f_divided_m%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        f_divided_m%fc(i,xyzt,s)= &
            a%fc(1,xyzt,s)*CONJG(b%mc(i,1))+ &
            a%fc(2,xyzt,s)*CONJG(b%mc(i,2))
      END DO
      END DO
      END DO
   END FUNCTION f_divided_m

!//
   FUNCTION g_predivided_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_predivided_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_predivided_g%parity=a%parity
      ELSE
        g_predivided_g%parity=-1
      ENDIF
      IF(a%dir.EQ.b%dir) THEN
        g_predivided_g%dir=a%dir
      ELSE
        g_predivided_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_predivided_g%fc(i,j,xyzt)= &
            CONJG(a%fc(1,i,xyzt))*b%fc(1,j,xyzt)+ &
            CONJG(a%fc(2,i,xyzt))*b%fc(2,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION g_predivided_g

   FUNCTION g_predivided_m(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(matrix), INTENT(IN) :: b
      TYPE(gauge_field) g_predivided_m
      INTEGER xyzt,i,j
      g_predivided_m%parity=a%parity
      g_predivided_m%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        g_predivided_m%fc(i,j,xyzt)= &
            CONJG(a%fc(1,i,xyzt))*b%mc(1,j)+ &
            CONJG(a%fc(2,i,xyzt))*b%mc(2,j)
      END DO
      END DO
      END DO
   END FUNCTION g_predivided_m

   FUNCTION m_predivided_g(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) m_predivided_g
      INTEGER xyzt,i,j
      m_predivided_g%parity=b%parity
      m_predivided_g%dir=b%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        m_predivided_g%fc(i,j,xyzt)= &
            CONJG(a%mc(1,i))*b%fc(1,j,xyzt)+ &
            CONJG(a%mc(2,i))*b%fc(2,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION m_predivided_g

   FUNCTION g_predivided_f(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) g_predivided_f
      INTEGER s,xyzt,i
      IF(a%parity.EQ.b%parity) THEN
        g_predivided_f%parity=a%parity
      ELSE
        g_predivided_f%parity=-1
      ENDIF
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        g_predivided_f%fc(i,xyzt,s)= &
            CONJG(a%fc(1,i,xyzt))*b%fc(1,xyzt,s)+ &
            CONJG(a%fc(2,i,xyzt))*b%fc(2,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION g_predivided_f

   FUNCTION m_predivided_f(a,b)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) m_predivided_f
      INTEGER s,xyzt,i
      m_predivided_f%parity=b%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        m_predivided_f%fc(i,xyzt,s)= &
            CONJG(a%mc(1,i))*b%fc(1,xyzt,s)+ &
            CONJG(a%mc(2,i))*b%fc(2,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION m_predivided_f

   FUNCTION f_diadic_f(a,b)
      TYPE(fermi_field), INTENT(IN) :: a,b
      TYPE(gauge_field) f_diadic_f
      INTEGER s,xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        f_diadic_f%parity=a%parity
      ELSE
        f_diadic_f%parity=-1
      ENDIF
      f_diadic_f%dir=0
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        f_diadic_f%fc(i,j,xyzt)= &
                 a%fc(i,xyzt,1)*CONJG(b%fc(j,xyzt,1))
        DO s=2,4
          f_diadic_f%fc(i,j,xyzt)=f_diadic_f%fc(i,j,xyzt)+ &
                   a%fc(i,xyzt,s)*CONJG(b%fc(j,xyzt,s)) 
        END DO
      END DO
      END DO
      END DO
   END FUNCTION f_diadic_f

!Dot
   FUNCTION g_dot_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(real_field) g_dot_g
      INTEGER xyzt,i,j
      IF(a%parity.EQ.b%parity) THEN
        g_dot_g%parity=a%parity
      ELSE
        g_dot_g%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
        g_dot_g%fc(xyzt)= &
          REAL(a%fc(1,1,xyzt),REAL8)*REAL(b%fc(1,1,xyzt),REAL8)+ &
          AIMAG(a%fc(1,1,xyzt))*AIMAG(b%fc(1,1,xyzt))+ &
          REAL(a%fc(2,1,xyzt),REAL8)*REAL(b%fc(2,1,xyzt),REAL8)+ &
          AIMAG(a%fc(2,1,xyzt))*AIMAG(b%fc(2,1,xyzt))
      DO j=2,2
        g_dot_g%fc(xyzt)=g_dot_g%fc(xyzt)+ &
          REAL(a%fc(1,j,xyzt),REAL8)*REAL(b%fc(1,j,xyzt),REAL8)+ &
          AIMAG(a%fc(1,j,xyzt))*AIMAG(b%fc(1,j,xyzt))+ &
          REAL(a%fc(2,j,xyzt),REAL8)*REAL(b%fc(2,j,xyzt),REAL8)+ &
          AIMAG(a%fc(2,j,xyzt))*AIMAG(b%fc(2,j,xyzt))
      END DO
      END DO
   END FUNCTION g_dot_g

!Gamma
   FUNCTION gamma_f(a,b)
      INTEGER, INTENT(IN) :: a  
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) :: gamma_f
      INTEGER :: xyzt,i
      gamma_f%parity=b%parity
      SELECT CASE(a) 
      CASE(1)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,4)=b%fc(i,xyzt,1)
        END DO
        END DO
      CASE(2)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,4)), &
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), & 
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE(3)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,2)=-b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,2)
        END DO
        END DO
      CASE(4)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,3)=-b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,4)
        END DO
        END DO
      CASE(5)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,4)), & 
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(11,22,33,44)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=(0._REAL8,0._REAL8)
          gamma_f%fc(i,xyzt,2)=(0._REAL8,0._REAL8)
          gamma_f%fc(i,xyzt,3)=(0._REAL8,0._REAL8)
          gamma_f%fc(i,xyzt,4)=(0._REAL8,0._REAL8)
        END DO
        END DO
      CASE(51)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,1)), & 
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), &
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), &
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(15)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), & 
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,4)), &
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(52)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=-b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,3)
        END DO
        END DO
      CASE(25)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,2)=-b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,3)=-b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,4)=b%fc(i,xyzt,3)
        END DO
        END DO
      CASE(53)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,1)), &
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), & 
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), &
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,4)), &
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
        END DO
        END DO
      CASE(35)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,2)), & 
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), &
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
        END DO
        END DO
      CASE(54)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), &
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), & 
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(45)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,4)), & 
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,1)), &
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(12)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=-b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,3)=-b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,4)=b%fc(i,xyzt,4)
        END DO
        END DO
      CASE(21)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,2)=-b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,4)
        END DO
        END DO
      CASE(23)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=-b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,2)=-b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,3)=-b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,3)
        END DO
        END DO
      CASE(32)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,1)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,4)=b%fc(i,xyzt,3)
        END DO
        END DO
      CASE(34)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), & 
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(43)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), &
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,4)), & 
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,1)), &
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(13)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), & 
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,4)), &
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), &
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(31)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,1)), & 
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), &
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,3)), &
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(24)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=-b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,2)=b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,3)=b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,4)=-b%fc(i,xyzt,1)
        END DO
        END DO
      CASE(42)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=b%fc(i,xyzt,4)
          gamma_f%fc(i,xyzt,2)=-b%fc(i,xyzt,3)
          gamma_f%fc(i,xyzt,3)=-b%fc(i,xyzt,2)
          gamma_f%fc(i,xyzt,4)=b%fc(i,xyzt,1)
        END DO
        END DO
      CASE(14)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX( AIMAG(b%fc(i,xyzt,4)), &
                                      -REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX( AIMAG(b%fc(i,xyzt,3)), & 
                                      -REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX(-AIMAG(b%fc(i,xyzt,2)), &
                                       REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX(-AIMAG(b%fc(i,xyzt,1)), &
                                       REAL(b%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE(41)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          gamma_f%fc(i,xyzt,1)=CMPLX(-AIMAG(b%fc(i,xyzt,4)), &
                                       REAL(b%fc(i,xyzt,4),REAL8),REAL8)
          gamma_f%fc(i,xyzt,2)=CMPLX(-AIMAG(b%fc(i,xyzt,3)), & 
                                       REAL(b%fc(i,xyzt,3),REAL8),REAL8)
          gamma_f%fc(i,xyzt,3)=CMPLX( AIMAG(b%fc(i,xyzt,2)), &
                                      -REAL(b%fc(i,xyzt,2),REAL8),REAL8)
          gamma_f%fc(i,xyzt,4)=CMPLX( AIMAG(b%fc(i,xyzt,1)), &
                                      -REAL(b%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE DEFAULT
        PRINT *,'index out of range for .Gamma. '
        STOP
      END SELECT
   END FUNCTION gamma_f

   FUNCTION f_gamma(a,b)
      TYPE(fermi_field), INTENT(IN) :: a
      INTEGER, INTENT(IN) :: b 
      TYPE(fermi_field) :: f_gamma
      INTEGER :: xyzt,i
      f_gamma%parity=a%parity
      SELECT CASE(b)
      CASE(1)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,4)=a%fc(i,xyzt,1)
        END DO
        END DO
      CASE(2)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), &
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,3)), & 
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE(3)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,2)=-a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,2)
        END DO
        END DO
      CASE(4)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,3)=-a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,4)
        END DO
        END DO
      CASE(5)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), & 
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(11,22,33,44)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=(0._REAL8,0._REAL8)
          f_gamma%fc(i,xyzt,2)=(0._REAL8,0._REAL8)
          f_gamma%fc(i,xyzt,3)=(0._REAL8,0._REAL8)
          f_gamma%fc(i,xyzt,4)=(0._REAL8,0._REAL8)
        END DO
        END DO
      CASE(51)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,1)), & 
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), &
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(15)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), & 
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,4)), &
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,3)), &
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(25)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=-a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,3)
        END DO
        END DO
      CASE(52)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,2)=-a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,3)=-a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,4)=a%fc(i,xyzt,3)
        END DO
        END DO
      CASE(53)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), & 
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,4)), &
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
        END DO
        END DO
      CASE(35)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), &
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,2)), & 
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,3)), &
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), &
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
        END DO
        END DO
      CASE(54)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), & 
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), &
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(45)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,3)), &
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,4)), & 
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(12)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=-a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,3)=-a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,4)=a%fc(i,xyzt,4)
        END DO
        END DO
      CASE(21)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,2)=-a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,4)
        END DO
        END DO
      CASE(23)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=-a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,2)=-a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,3)=-a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,3)
        END DO
        END DO
      CASE(32)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,1)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,4)=a%fc(i,xyzt,3)
        END DO
        END DO
      CASE(43)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,3)), &
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), & 
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), &
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(34)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,4)), & 
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
        END DO
        END DO
      CASE(31)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), & 
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,4)), &
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), &
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(13)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,1)), & 
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), &
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,3)), &
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
        END DO
        END DO
      CASE(24)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=-a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,2)=a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,3)=a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,4)=-a%fc(i,xyzt,1)
        END DO
        END DO
      CASE(42)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=a%fc(i,xyzt,4)
          f_gamma%fc(i,xyzt,2)=-a%fc(i,xyzt,3)
          f_gamma%fc(i,xyzt,3)=-a%fc(i,xyzt,2)
          f_gamma%fc(i,xyzt,4)=a%fc(i,xyzt,1)
        END DO
        END DO
      CASE(41)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX( AIMAG(a%fc(i,xyzt,4)), &
                                      -REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX( AIMAG(a%fc(i,xyzt,3)), & 
                                      -REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX(-AIMAG(a%fc(i,xyzt,2)), &
                                       REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX(-AIMAG(a%fc(i,xyzt,1)), &
                                       REAL(a%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE(14)
        DO xyzt=0,NXYZT2-1
        DO i=1,2
          f_gamma%fc(i,xyzt,1)=CMPLX(-AIMAG(a%fc(i,xyzt,4)), &
                                       REAL(a%fc(i,xyzt,4),REAL8),REAL8)
          f_gamma%fc(i,xyzt,2)=CMPLX(-AIMAG(a%fc(i,xyzt,3)), & 
                                       REAL(a%fc(i,xyzt,3),REAL8),REAL8)
          f_gamma%fc(i,xyzt,3)=CMPLX( AIMAG(a%fc(i,xyzt,2)), &
                                      -REAL(a%fc(i,xyzt,2),REAL8),REAL8)
          f_gamma%fc(i,xyzt,4)=CMPLX( AIMAG(a%fc(i,xyzt,1)), &
                                      -REAL(a%fc(i,xyzt,1),REAL8),REAL8)
        END DO
        END DO
      CASE DEFAULT
        PRINT *,'index out of range for .Gamma. '
        STOP
      END SELECT
   END FUNCTION f_gamma

!Lambda
   FUNCTION lambda_g(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) :: lambda_g
      INTEGER :: xyzt
      lambda_g%parity=b%parity
      lambda_g%dir=b%dir
      SELECT CASE(a) 
                          !a=1 is \lambda_1*g, etc...
      CASE(1)
        DO xyzt=0,NXYZT2-1
          lambda_g%fc(1,1,xyzt)=b%fc(2,1,xyzt)
          lambda_g%fc(1,2,xyzt)=b%fc(2,2,xyzt)
          lambda_g%fc(2,1,xyzt)=b%fc(1,1,xyzt)
          lambda_g%fc(2,2,xyzt)=b%fc(1,2,xyzt)
        END DO
      CASE(2)
        DO xyzt=0,NXYZT2-1
          lambda_g%fc(1,1,xyzt)=CMPLX(AIMAG(b%fc(2,1,xyzt)),&
                         -REAL(b%fc(2,1,xyzt),REAL8),REAL8)
          lambda_g%fc(1,2,xyzt)=CMPLX(AIMAG(b%fc(2,2,xyzt)),&
                         -REAL(b%fc(2,2,xyzt),REAL8),REAL8)
          lambda_g%fc(2,1,xyzt)=CMPLX(-AIMAG(b%fc(1,1,xyzt)),&
                         REAL(b%fc(1,1,xyzt),REAL8),REAL8)
          lambda_g%fc(2,2,xyzt)=CMPLX(-AIMAG(b%fc(1,2,xyzt)),&
                         REAL(b%fc(1,2,xyzt),REAL8),REAL8)
        END DO
      CASE(3)
        DO xyzt=0,NXYZT2-1
          lambda_g%fc(1,1,xyzt)=b%fc(1,1,xyzt)
          lambda_g%fc(1,2,xyzt)=b%fc(1,2,xyzt)
          lambda_g%fc(2,1,xyzt)=-b%fc(2,1,xyzt)
          lambda_g%fc(2,2,xyzt)=-b%fc(2,2,xyzt)
        END DO
      END SELECT
   END FUNCTION lambda_g

   FUNCTION g_lambda(a,b)
      TYPE(gauge_field), INTENT(IN) :: a
      INTEGER, INTENT(IN) :: b
      TYPE(gauge_field) :: g_lambda
      INTEGER :: xyzt
      g_lambda%parity=a%parity
      g_lambda%dir=a%dir
      SELECT CASE(b) 
                          !b=1 is \g*lambda_1, etc...
      CASE(1)
        DO xyzt=0,NXYZT2-1
          g_lambda%fc(1,1,xyzt)=a%fc(1,2,xyzt)
          g_lambda%fc(2,1,xyzt)=a%fc(2,2,xyzt)
          g_lambda%fc(1,2,xyzt)=a%fc(1,1,xyzt)
          g_lambda%fc(2,2,xyzt)=a%fc(2,1,xyzt)
        END DO
      CASE(2)
        DO xyzt=0,NXYZT2-1
          g_lambda%fc(1,1,xyzt)=CMPLX(-AIMAG(a%fc(1,2,xyzt)),&
                          REAL(a%fc(1,2,xyzt),REAL8),REAL8)
          g_lambda%fc(2,1,xyzt)=CMPLX(-AIMAG(a%fc(2,2,xyzt)),&
                          REAL(a%fc(2,2,xyzt),REAL8),REAL8)
          g_lambda%fc(1,2,xyzt)=CMPLX(AIMAG(a%fc(1,1,xyzt)),&
                         -REAL(a%fc(1,1,xyzt),REAL8),REAL8)
          g_lambda%fc(2,2,xyzt)=CMPLX(AIMAG(a%fc(2,1,xyzt)),&
                         -REAL(a%fc(2,1,xyzt),REAL8),REAL8)
        END DO
      CASE(3)
        DO xyzt=0,NXYZT2-1
          g_lambda%fc(1,1,xyzt)=a%fc(1,1,xyzt)
          g_lambda%fc(2,1,xyzt)=a%fc(2,1,xyzt)
          g_lambda%fc(1,2,xyzt)=-a%fc(1,2,xyzt)
          g_lambda%fc(2,2,xyzt)=-a%fc(2,2,xyzt)
        END DO
      END SELECT
   END FUNCTION g_lambda


!I
   FUNCTION I_times_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(gauge_field) :: I_times_g
      INTEGER xyzt,i,j
      I_times_g%parity=a%parity
      I_times_g%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        I_times_g%fc(i,j,xyzt)= &
             CMPLX(-AIMAG(a%fc(i,j,xyzt)), &
                     REAL(a%fc(i,j,xyzt),REAL8),REAL8)
      END DO
      END DO
      END DO
   END FUNCTION I_times_g

   FUNCTION I_times_f(a)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(fermi_field) :: I_times_f
      INTEGER xyzt,i,s
      I_times_f%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        I_times_f%fc(i,xyzt,s)= &
             CMPLX(-AIMAG(a%fc(i,xyzt,s)), &
                     REAL(a%fc(i,xyzt,s),REAL8),REAL8) 
      END DO
      END DO
      END DO
   END FUNCTION I_times_f

   FUNCTION I_times_c(a)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(complex_field) :: I_times_c
      INTEGER xyzt
      I_times_c%parity=a%parity
      DO xyzt=0,NXYZT2-1
        I_times_c%fc(xyzt)= &
             CMPLX(-AIMAG(a%fc(xyzt)), &
                     REAL(a%fc(xyzt),REAL8),REAL8)
      END DO
   END FUNCTION I_times_c

   FUNCTION I_times_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(complex_field) :: I_times_r
      INTEGER xyzt
      I_times_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        I_times_r%fc(xyzt)=CMPLX(0.0_REAL8,a%fc(xyzt),REAL8)
      END DO
   END FUNCTION I_times_r


!Minus
   FUNCTION minus_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(gauge_field) minus_g
      INTEGER xyzt,i,j
      minus_g%parity=a%parity
      minus_g%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        minus_g%fc(i,j,xyzt)=-a%fc(i,j,xyzt)
      END DO
      END DO
      END DO
   END FUNCTION minus_g

   FUNCTION minus_f(a)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(fermi_field) minus_f
      INTEGER s,xyzt,i
      minus_f%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        minus_f%fc(i,xyzt,s)=-a%fc(i,xyzt,s)
      END DO
      END DO
      END DO
   END FUNCTION minus_f

   FUNCTION minus_c(a)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(complex_field) minus_c
      INTEGER xyzt
      minus_c%parity=a%parity
      DO xyzt=0,NXYZT2-1
        minus_c%fc(xyzt)=-a%fc(xyzt)
      END DO
   END FUNCTION minus_c

   FUNCTION minus_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) minus_r
      INTEGER xyzt
      minus_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        minus_r%fc(xyzt)=-a%fc(xyzt)
      END DO
   END FUNCTION minus_r

   FUNCTION minus_ge(a)
      TYPE(generator_field), INTENT(IN) :: a
      TYPE(generator_field) minus_ge
      INTEGER xyzt,i
      minus_ge%parity=a%parity
      minus_ge%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO i=1,3
        minus_ge%fc(i,xyzt)=-a%fc(i,xyzt)
      END DO
      END DO
   END FUNCTION minus_ge

!Conjg
   FUNCTION conjg_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(gauge_field) conjg_g
      INTEGER xyzt,i,j
      conjg_g%parity=a%parity
      conjg_g%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        conjg_g%fc(i,j,xyzt)=CONJG(a%fc(i,j,xyzt))
      END DO
      END DO
      END DO
   END FUNCTION conjg_g

   FUNCTION conjg_f(a)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(fermi_field) conjg_f
      INTEGER s,xyzt,i
      conjg_f%parity=a%parity
      DO s=1,4
      DO xyzt=0,NXYZT2-1
      DO i=1,2
        conjg_f%fc(i,xyzt,s)=CONJG(a%fc(i,xyzt,s))
      END DO
      END DO
      END DO
   END FUNCTION conjg_f

   FUNCTION conjg_c(a)
      TYPE(complex_field), INTENT(IN) :: a
      TYPE(complex_field) conjg_c
      INTEGER xyzt
      conjg_c%parity=a%parity
      DO xyzt=0,NXYZT2-1
        conjg_c%fc(xyzt)=CONJG(a%fc(xyzt))
      END DO
   END FUNCTION conjg_c

   FUNCTION conjg_m(a)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(matrix) conjg_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        conjg_m%mc(i,j)=CONJG(a%mc(i,j))
      END DO
      END DO
   END FUNCTION conjg_m

!Adj
   FUNCTION adj_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(gauge_field) adj_g
      INTEGER xyzt,i,j
      adj_g%parity=a%parity
      adj_g%dir=a%dir
      DO xyzt=0,NXYZT2-1
      DO j=1,2
      DO i=1,2
        adj_g%fc(i,j,xyzt)=CONJG(a%fc(j,i,xyzt))
      END DO
      END DO
      END DO
   END FUNCTION adj_g

   FUNCTION adj_m(a)
      TYPE(matrix), INTENT(IN) :: a
      TYPE(matrix) adj_m
      INTEGER i,j
      DO j=1,2
      DO i=1,2
        adj_m%mc(i,j)=CONJG(a%mc(j,i))
      END DO
      END DO
   END FUNCTION adj_m

!Ctr
   FUNCTION ctr_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(complex_field) ctr_g
      INTEGER xyzt
      ctr_g%parity=a%parity
      DO xyzt=0,NXYZT2-1
        ctr_g%fc(xyzt)= a%fc(1,1,xyzt)+ &
                        a%fc(2,2,xyzt)
      END DO
   END FUNCTION ctr_g

   FUNCTION ctr_m(a)
      TYPE(matrix), INTENT(IN) :: a
      COMPLEX(REAL8) ctr_m
      ctr_m=a%mc(1,1)+a%mc(2,2)
   END FUNCTION ctr_m

!Tr
   FUNCTION tr_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(real_field) tr_g
      INTEGER xyzt
      tr_g%parity=a%parity
      DO xyzt=0,NXYZT2-1
        tr_g%fc(xyzt)=REAL(a%fc(1,1,xyzt),REAL8)+ &
                      REAL(a%fc(2,2,xyzt),REAL8)
      END DO
   END FUNCTION tr_g

   FUNCTION tr_m(a)
      TYPE(matrix), INTENT(IN) :: a
      REAL(REAL8) tr_m
      tr_m=REAL(a%mc(1,1)+a%mc(2,2),REAL8)
   END FUNCTION tr_m
 
!Sqrt
   FUNCTION sqrt_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) sqrt_r
      INTEGER xyzt
      sqrt_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        context(xyzt)=a%fc(xyzt)>=0
        sqrt_r%fc(xyzt)=SQRT(ABS(a%fc(xyzt)))
      END DO
   END FUNCTION sqrt_r

!Exp
   FUNCTION exp_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) exp_r
      INTEGER xyzt
      exp_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        exp_r%fc(xyzt)=EXP(a%fc(xyzt))
      END DO
   END FUNCTION exp_r

!Log
   FUNCTION log_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) log_r
      INTEGER xyzt
      log_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        log_r%fc(xyzt)=LOG(a%fc(xyzt))
      END DO
   END FUNCTION log_r

!Cos
   FUNCTION cos_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) cos_r
      INTEGER xyzt
      cos_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        cos_r%fc(xyzt)=COS(a%fc(xyzt))
      END DO
   END FUNCTION cos_r

!Sin
   FUNCTION sin_r(a)
      TYPE(real_field), INTENT(IN) :: a
      TYPE(real_field) sin_r
      INTEGER xyzt
      sin_r%parity=a%parity
      DO xyzt=0,NXYZT2-1
        sin_r%fc(xyzt)=SIN(a%fc(xyzt))
      END DO
   END FUNCTION sin_r


END MODULE

