
!  Program qcdf90, module assign_mixed, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.
 
! su2 version  2006.08.28  by y.k.

MODULE assign_mixed

USE precisions
USE constants
USE global_module

   IMPLICIT NONE  


   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE g_asgn_u,g_asgn_ge,g_asgn_m, &
                       g_asgn_complex,g_asgn_real, &
                       u_asgn_g,f_asgn_g,f_asgn_c,f_asgn_r, &
                       f_asgn_complex, f_asgn_real, &
                       c_asgn_g,c_asgn_r, &
                       c_asgn_complex,c_asgn_real,r_asgn_g,r_asgn_f, &
                       r_asgn_real,ge_asgn_g, &
                       ge_asgn_r,ge_asgn_real, &
                       complex_asgn_c,real_asgn_c,real_asgn_r
   END INTERFACE

CONTAINS

   SUBROUTINE g_asgn_u(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      TYPE(full_gauge_field), INTENT(IN) :: b
      INTEGER i,j,xyzt
      SELECT CASE(assign_type)

      CASE('=')
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=b%uc(a%parity,a%dir)%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('t')
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=b%uc(a%parity,a%dir)%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE('f')
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=b%uc(a%parity,a%dir)%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_u

   SUBROUTINE g_asgn_ge(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      TYPE(generator_field), INTENT(IN) :: b
      REAL(REAL8) s,c1,c2
      COMPLEX(REAL8) ck1,ck2,ck3,ck4
      INTEGER xyzt,i,j,k
      COMPLEX(REAL8), DIMENSION(2,2) :: ms
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         a%dir=b%dir            
         DO xyzt=0,NXYZT2-1
           a%fc(1,1,xyzt)= b%fc(3,xyzt)
           a%fc(1,2,xyzt)=CMPLX(b%fc(1,xyzt),-b%fc(2,xyzt),REAL8)
           a%fc(2,1,xyzt)=CMPLX(b%fc(1,xyzt), b%fc(2,xyzt),REAL8)
           a%fc(2,2,xyzt)=-b%fc(3,xyzt)
         END DO

      CASE('E')
        a%parity=b%parity    
        a%dir=b%dir            
        DO xyzt=0,NXYZT2-1

           s=b%fc(1,xyzt)**2
           DO i=2,3
              s=s+b%fc(i,xyzt)**2
           END DO
           s=SQRT(s)

           c1=COS(s)
           IF(ABS(s)>0.00000001_REAL8) THEN
              c2=SIN(s)/s
           ELSE
              c2=1._REAL8
           ENDIF 

!  ms=c1*UNIT+IU*c2*(.Matrix.hs.), inlined:
     

         a%fc(1,1,xyzt)=c1+c2*CMPLX(0._REAL8,  b%fc(3,xyzt),REAL8)
         a%fc(2,2,xyzt)=c1+c2*CMPLX(0._REAL8, -b%fc(3,xyzt),REAL8)
         a%fc(1,2,xyzt)=c2*CMPLX( b%fc(2,xyzt),b%fc(1,xyzt),REAL8)
         a%fc(2,1,xyzt)=c2*CMPLX(-b%fc(2,xyzt),b%fc(1,xyzt),REAL8)

        END DO         

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_ge

   SUBROUTINE g_asgn_m(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      TYPE(matrix), INTENT(IN) :: b
      INTEGER i,j,xyzt
      COMPLEX(REAL8) aux(2,2)
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=-1
         a%dir=0
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=b%mc(i,j)
         END DO
         END DO
         END DO

      CASE('*')
         SELECT CASE(assign_spec)

         CASE(0)
            DO xyzt=0,NXYZT2-1
               DO j=1,2
               DO i=1,2
                  aux(i,j)=a%fc(i,1,xyzt)*b%mc(1,j)+ &
                           a%fc(i,2,xyzt)*b%mc(2,j)
               END DO
               END DO
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=aux(i,j)
               END DO
               END DO
            END DO

         CASE(-1)
            DO xyzt=0,NXYZT2-1
               DO j=1,2
               DO i=1,2
                  aux(i,j)=b%mc(i,1)*a%fc(1,j,xyzt)+ &
                           b%mc(i,2)*a%fc(2,j,xyzt)
               END DO
               END DO
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=aux(i,j)
               END DO
               END DO
            END DO

         CASE DEFAULT
           PRINT '("Illegal value for assign_spec")'
           STOP
         END SELECT

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_m

   SUBROUTINE g_asgn_complex(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      INTEGER i,j,xyzt
      SELECT CASE(assign_type)

      CASE('*')
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
             a%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b      
         END DO
         END DO
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
             a%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b      
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_complex

   SUBROUTINE g_asgn_real(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      REAL(REAL8), INTENT(IN) :: b
      INTEGER i,j,xyzt
      SELECT CASE(assign_type)

      CASE('*')
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
             a%fc(i,j,xyzt)=a%fc(i,j,xyzt)*b      
         END DO
         END DO
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
             a%fc(i,j,xyzt)=a%fc(i,j,xyzt)/b      
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_real

   SUBROUTINE u_asgn_g(a,b)
      TYPE(full_gauge_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      INTEGER i,j,xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%uc(b%parity,b%dir)%parity=b%parity
         a%uc(b%parity,b%dir)%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%uc(b%parity,b%dir)%fc(i,j,xyzt)=b%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('t')
         a%uc(b%parity,b%dir)%parity=b%parity
         a%uc(b%parity,b%dir)%dir=b%dir
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%uc(b%parity,b%dir)%fc(i,j,xyzt)=b%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE('f')
         a%uc(b%parity,b%dir)%parity=b%parity
         a%uc(b%parity,b%dir)%dir=b%dir
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%uc(b%parity,b%dir)%fc(i,j,xyzt)=b%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE u_asgn_g

   SUBROUTINE f_asgn_g(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      INTEGER i,xyzt,s
      COMPLEX(REAL8) auxil(2)
      SELECT CASE(assign_type)

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         SELECT CASE(assign_spec)

         CASE(0)
            DO s=1,4
            DO xyzt=0,NXYZT2-1
               DO i=1,2
                  auxil(i)=a%fc(1,xyzt,s)*b%fc(i,1,xyzt)+ &
                           a%fc(2,xyzt,s)*b%fc(i,2,xyzt)
               END DO
               DO i=1,2
                  a%fc(i,xyzt,s)=auxil(i)
               END DO
            END DO
            END DO

         CASE(-1)
            DO s=1,4
            DO xyzt=0,NXYZT2-1
               DO i=1,2
                  auxil(i)=a%fc(1,xyzt,s)*b%fc(1,i,xyzt)+ &
                           a%fc(2,xyzt,s)*b%fc(2,i,xyzt)
               END DO
               DO i=1,2
                  a%fc(i,xyzt,s)=auxil(i)
               END DO
            END DO
            END DO

         CASE DEFAULT
           PRINT '("Illegal value for assign_spec")'
           STOP
         END SELECT

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1

         SELECT CASE(assign_spec)

         CASE(0)
            DO s=1,4
            DO xyzt=0,NXYZT2-1
               DO i=1,2
                  auxil(i)=a%fc(1,xyzt,s)*CONJG(b%fc(i,1,xyzt))+ &
                           a%fc(2,xyzt,s)*CONJG(b%fc(i,2,xyzt))
               END DO
               DO i=1,2
                  a%fc(i,xyzt,s)=auxil(i)
               END DO
            END DO
            END DO

         CASE(-1)
            DO s=1,4
            DO xyzt=0,NXYZT2-1
               DO i=1,2
                  auxil(i)=CONJG(b%fc(1,i,xyzt))*a%fc(1,xyzt,s)+ &
                           CONJG(b%fc(2,i,xyzt))*a%fc(2,xyzt,s)
               END DO
               DO i=1,2
                  a%fc(i,xyzt,s)=auxil(i)
               END DO
            END DO
            END DO

         CASE DEFAULT
           PRINT '("Illegal value for assign_spec")'
           STOP
         END SELECT

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_g

   SUBROUTINE f_asgn_c(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      TYPE(complex_field), INTENT(IN) :: b
      INTEGER i,xyzt,s
      SELECT CASE(assign_type)

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
         END DO
         END DO
         END DO

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b%fc(xyzt)
!???        a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_c

   SUBROUTINE f_asgn_r(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      TYPE(real_field), INTENT(IN) :: b
      INTEGER i,xyzt,s
      SELECT CASE(assign_type)

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
         END DO
         END DO
         END DO

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b%fc(xyzt)
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_r

   SUBROUTINE f_asgn_complex(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      INTEGER i,xyzt,s
      SELECT CASE(assign_type)

      CASE('*')
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b
         END DO
         END DO
         END DO

      CASE('/')
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_complex


   SUBROUTINE f_asgn_real(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      REAL(REAL8), INTENT(IN) :: b
      INTEGER i,xyzt,s
      SELECT CASE(assign_type)

      CASE('*')
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)*b
         END DO
         END DO
         END DO

      CASE('/')
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=a%fc(i,xyzt,s)/b
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_real

   SUBROUTINE c_asgn_g(a,b)
      TYPE(complex_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
           a%fc(xyzt)= b%fc(1,1,xyzt)+b%fc(2,2,xyzt)
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE c_asgn_g

   SUBROUTINE c_asgn_r(a,b)
      TYPE(complex_field), INTENT(INOUT) :: a
      TYPE(real_field), INTENT(IN) :: b
      INTEGER xyzt
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

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE c_asgn_r


   SUBROUTINE c_asgn_complex(a,b)
      TYPE(complex_field), INTENT(INOUT) :: a
      COMPLEX(REAL8), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=b
         END DO

      CASE('+')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)+b
         END DO

      CASE('-')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)-b
         END DO

      CASE('*')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)*b
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)/b
         END DO

      CASE('M')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=-b
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE c_asgn_complex

   SUBROUTINE c_asgn_real(a,b)
      TYPE(complex_field), INTENT(INOUT) :: a
      REAL(REAL8), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=b
         END DO

      CASE('+')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)+b
         END DO

      CASE('-')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)-b
         END DO

      CASE('*')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)*b
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)/b
         END DO

      CASE('M')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=-b
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE c_asgn_real

   SUBROUTINE r_asgn_f(a,b)
      TYPE(real_field), INTENT(INOUT) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
         a%fc(xyzt)= &
            REAL(b%fc(1,xyzt,1),REAL8)*REAL(b%fc(1,xyzt,1),REAL8)+&
           AIMAG(b%fc(1,xyzt,1))*AIMAG(b%fc(1,xyzt,1))+&
            REAL(b%fc(1,xyzt,2),REAL8)*REAL(b%fc(1,xyzt,2),REAL8)+&
           AIMAG(b%fc(1,xyzt,2))*AIMAG(b%fc(1,xyzt,2))+&
            REAL(b%fc(1,xyzt,3),REAL8)*REAL(b%fc(1,xyzt,3),REAL8)+&
           AIMAG(b%fc(1,xyzt,3))*AIMAG(b%fc(1,xyzt,3))+&
            REAL(b%fc(1,xyzt,4),REAL8)*REAL(b%fc(1,xyzt,4),REAL8)+&
           AIMAG(b%fc(1,xyzt,4))*AIMAG(b%fc(1,xyzt,4))+&
            REAL(b%fc(2,xyzt,1),REAL8)*REAL(b%fc(2,xyzt,1),REAL8)+&
           AIMAG(b%fc(2,xyzt,1))*AIMAG(b%fc(2,xyzt,1))+&
            REAL(b%fc(2,xyzt,2),REAL8)*REAL(b%fc(2,xyzt,2),REAL8)+&
           AIMAG(b%fc(2,xyzt,2))*AIMAG(b%fc(2,xyzt,2))+&
            REAL(b%fc(2,xyzt,3),REAL8)*REAL(b%fc(2,xyzt,3),REAL8)+&
           AIMAG(b%fc(2,xyzt,3))*AIMAG(b%fc(2,xyzt,3))+&
            REAL(b%fc(2,xyzt,4),REAL8)*REAL(b%fc(2,xyzt,4),REAL8)+&
           AIMAG(b%fc(2,xyzt,4))*AIMAG(b%fc(2,xyzt,4))
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE r_asgn_f

   SUBROUTINE r_asgn_g(a,b)
      TYPE(real_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1        
            a%fc(xyzt)=REAL(b%fc(1,1,xyzt),REAL8)+ &
                       REAL(b%fc(2,2,xyzt),REAL8)
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE r_asgn_g

   SUBROUTINE r_asgn_c(a,b)
      TYPE(real_field), INTENT(INOUT) :: a
      TYPE(complex_field), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=REAL(b%fc(xyzt),REAL8)
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE r_asgn_c

   SUBROUTINE r_asgn_real(a,b)
      TYPE(real_field), INTENT(INOUT) :: a
      REAL(REAL8), INTENT(IN) :: b
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
             a%fc(xyzt)=b
         END DO

      CASE('+')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)+b
         END DO

      CASE('-')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)-b
         END DO

      CASE('*')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)*b
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=a%fc(xyzt)/b
         END DO

      CASE('M')
         a%parity=-1
         DO xyzt=0,NXYZT2-1
            a%fc(xyzt)=-b
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE r_asgn_real

   SUBROUTINE ge_asgn_g(a,b)
      TYPE(generator_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      REAL(REAL8) aux_sq(3)
      INTEGER xyzt
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
           a%fc(1,xyzt)=0.5_REAL8*AIMAG(b%fc(1,2,xyzt)+b%fc(2,1,xyzt))
           a%fc(2,xyzt)=0.5_REAL8*REAL(b%fc(1,2,xyzt)-b%fc(2,1,xyzt),REAL8)
           a%fc(3,xyzt)=0.5_REAL8*AIMAG(b%fc(1,1,xyzt)-b%fc(2,2,xyzt))
      END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE ge_asgn_g

   SUBROUTINE ge_asgn_r(a,b)
      TYPE(generator_field), INTENT(INOUT) :: a
      TYPE(real_field), INTENT(IN) :: b
      INTEGER i,xyzt
      SELECT CASE(assign_type)

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)*b%fc(xyzt)
         END DO
         END DO

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)/b%fc(xyzt)
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE ge_asgn_r

   SUBROUTINE ge_asgn_real(a,b)
      TYPE(generator_field), INTENT(INOUT) :: a
      REAL(REAL8), INTENT(IN) :: b
      INTEGER i,xyzt
      SELECT CASE(assign_type)

      CASE('*')
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)*b
         END DO
         END DO

      CASE('/')
         DO xyzt=0,NXYZT2-1
         DO i=1,3
            a%fc(i,xyzt)=a%fc(i,xyzt)/b
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE ge_asgn_real

   SUBROUTINE complex_asgn_c(a,b)
      COMPLEX(REAL8), INTENT(INOUT) :: a
      TYPE(complex_field), INTENT(IN) :: b
      INTEGER xyzt
      a=(0._REAL8,0._REAL8)
      SELECT CASE(assign_type)

      CASE('=')
         DO xyzt=0,NXYZT2-1
            a=a+b%fc(xyzt)
         END DO

      CASE('t')
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               a=a+b%fc(xyzt)
            ENDIF
         END DO

      CASE('f')
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               a=a+b%fc(xyzt)
            ENDIF
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE complex_asgn_c

   SUBROUTINE real_asgn_c(a,b)
      REAL(REAL8), INTENT(INOUT) :: a
      TYPE(complex_field), INTENT(IN) :: b
      INTEGER xyzt
      a=0._REAL8
      SELECT CASE(assign_type)

      CASE('=')
         DO xyzt=0,NXYZT2-1
            a=a+REAL(b%fc(xyzt),REAL8)
         END DO

      CASE('t')
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               a=a+REAL(b%fc(xyzt),REAL8)
            ENDIF
         END DO

      CASE('f')
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               a=a+REAL(b%fc(xyzt),REAL8)
            ENDIF
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE real_asgn_c

   SUBROUTINE real_asgn_r(a,b)
      REAL(REAL8), INTENT(INOUT) :: a
      TYPE(real_field), INTENT(IN) :: b
      INTEGER xyzt
      a=0._REAL8
      SELECT CASE(assign_type)

      CASE('=')
         DO xyzt=0,NXYZT2-1
            a=a+b%fc(xyzt)
         END DO
 
      CASE('t')
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               a=a+b%fc(xyzt)
            ENDIF
         END DO

      CASE('f')
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               a=a+b%fc(xyzt)
            ENDIF
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE real_asgn_r

END MODULE

