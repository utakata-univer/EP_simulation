!  Program qcdf90, module assign_isotype3, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.
 
! su2 version  2006.08.28 by y.k.

MODULE assign_isotype3

USE precisions
USE constants
USE global_module

   IMPLICIT NONE  

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE u_asgn_u,c_asgn_c,r_asgn_r,ge_asgn_ge
   END INTERFACE

CONTAINS

   SUBROUTINE u_asgn_u(a,b)
      TYPE(full_gauge_field), INTENT(INOUT) :: a
      TYPE(full_gauge_field), INTENT(IN) :: b
      INTEGER i,j,xyzt,eo,m
      SELECT CASE(assign_type)

      CASE('=')
         DO m=1,4
         DO eo=0,1
            a%uc(eo,m)%parity=eo
            a%uc(eo,m)%dir=m
            DO xyzt=0,NXYZT2-1
            DO j=1,2
            DO i=1,2
               a%uc(eo,m)%fc(i,j,xyzt)=b%uc(eo,m)%fc(i,j,xyzt)
            END DO
            END DO
            END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE u_asgn_u

   SUBROUTINE c_asgn_c(a,b)
      TYPE(complex_field), INTENT(INOUT) :: a
      TYPE(complex_field), INTENT(IN) :: b
      INTEGER i,j,xyzt,k,m,eo
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=b%fc(xyzt)
         END DO

      CASE('+')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
         END DO

      CASE('-')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
         END DO

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
         END DO

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
         END DO

      CASE('C')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=CONJG(b%fc(xyzt))
         END DO

      CASE('I')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           a%fc(xyzt)=CMPLX(-AIMAG(b%fc(xyzt)), &
                              REAL(b%fc(xyzt),REAL8),REAL8) 
         END DO

      CASE('M')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           a%fc(xyzt)=-b%fc(xyzt)
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE c_asgn_c

   SUBROUTINE r_asgn_r(a,b)
      TYPE(real_field), INTENT(INOUT) :: a
      TYPE(real_field), INTENT(IN) :: b
      INTEGER i,j,xyzt,k,m,eo
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=b%fc(xyzt)
         END DO

      CASE('+')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)+b%fc(xyzt)
         END DO

      CASE('-')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)-b%fc(xyzt)
         END DO

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)*b%fc(xyzt)
         END DO

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)/b%fc(xyzt)
         END DO

      CASE('M')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           a%fc(xyzt)=-b%fc(xyzt)
         END DO

      CASE('R')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           context(xyzt)=b%fc(xyzt)>=0
           a%fc(xyzt)=SQRT(ABS(b%fc(xyzt)))
         END DO

      CASE('E')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           a%fc(xyzt)=EXP(b%fc(xyzt))
         END DO


      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE r_asgn_r

   SUBROUTINE ge_asgn_ge(a,b)
      TYPE(generator_field), INTENT(INOUT) :: a
      TYPE(generator_field), INTENT(IN) :: b
      REAL(REAL8) aux8,aux12,aux45,aux67
      REAL(REAL8) aux_sq(3)
      INTEGER i,xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=b%fc(i,xyzt)
         END DO
         END DO

      CASE('+')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)+b%fc(i,xyzt)
         END DO
         END DO

      CASE('-')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)-b%fc(i,xyzt)
         END DO
         END DO

      CASE('M')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=-b%fc(i,xyzt)
         END DO
         END DO

      CASE('S')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
            a%fc(1,xyzt)=0.0_REAL8
            a%fc(2,xyzt)=0.0_REAL8
            a%fc(3,xyzt)=0.0_REAL8
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE ge_asgn_ge

END MODULE
