
!  Program qcdf90, module conditionals, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

MODULE conditionals

USE precisions
USE global_module

   IMPLICIT NONE
      
                    ! The relational operators between real fields set
                    ! the global variable context to .TRUE. where the
                    ! relation is satisfied, to .FALSE. otherwise.
                    ! They return .TRUE. if the fields have the same
                    ! parity, .FALSE. if they have different parity
                    ! or their parity is undefined.

   INTERFACE OPERATOR(>)
      MODULE PROCEDURE r_greater_r, r_greater_real, real_greater_r
   END INTERFACE

   INTERFACE OPERATOR(>=)
      MODULE PROCEDURE r_greater_eq_r, r_greater_eq_real, real_greater_eq_r
   END INTERFACE

   INTERFACE OPERATOR(<)
      MODULE PROCEDURE r_less_r, r_less_real, real_less_r
   END INTERFACE

   INTERFACE OPERATOR(<=)
      MODULE PROCEDURE r_less_eq_r, r_less_eq_real, real_less_eq_r
   END INTERFACE

   INTERFACE OPERATOR(==)
      MODULE PROCEDURE r_equal_r, r_equal_real, real_equal_r
   END INTERFACE

   INTERFACE OPERATOR(/=)
      MODULE PROCEDURE r_not_equal_r, r_not_equal_real, real_not_equal_r
   END INTERFACE

                    ! .Xor., taken between any type of fields, returns
                    ! a field of the same type having as elements the
                    ! corresponding elements of the first operand where 
                    ! the global variable context is .TRUE., the elements 
                    ! of the second operand where context is .FALSE.
                    ! The parity of the returned field is the common
                    ! parity of the two operands if they have the same
                    ! parity, otherwise it is undefined. (For the gauge
                    ! fields the variable direction of the returned field 
                    ! is the common direction of the operands if they have 
                    ! the same direction, otherwise it is set to 0.    

   INTERFACE OPERATOR(.Xor.)  
      MODULE PROCEDURE r_xor_r, c_xor_c, f_xor_f, g_xor_g, ge_xor_ge
   END INTERFACE

CONTAINS

   FUNCTION r_greater_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_greater_r
      INTEGER xyzt
      r_greater_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)>b%fc(xyzt)
      END DO
   END FUNCTION r_greater_r   

   FUNCTION r_greater_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_greater_real
      INTEGER xyzt
      r_greater_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)>b
      END DO
   END FUNCTION r_greater_real   

   FUNCTION real_greater_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_greater_r
      INTEGER xyzt
      real_greater_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
           context(xyzt)=a>b%fc(xyzt)
      END DO
   END FUNCTION real_greater_r   

   FUNCTION r_greater_eq_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_greater_eq_r
      INTEGER xyzt
      r_greater_eq_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)>=b%fc(xyzt)
      END DO
   END FUNCTION r_greater_eq_r   

   FUNCTION r_greater_eq_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_greater_eq_real
      INTEGER xyzt
      r_greater_eq_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)>=b
      END DO
   END FUNCTION r_greater_eq_real   

   FUNCTION real_greater_eq_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_greater_eq_r
      INTEGER xyzt
      real_greater_eq_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a>=b%fc(xyzt)
      END DO
   END FUNCTION real_greater_eq_r   

   FUNCTION r_less_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_less_r
      INTEGER xyzt
      r_less_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)<b%fc(xyzt)
      END DO
   END FUNCTION r_less_r   

   FUNCTION r_less_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_less_real
      INTEGER xyzt
      r_less_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)<b
      END DO
   END FUNCTION r_less_real   

   FUNCTION real_less_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_less_r
      INTEGER xyzt
      real_less_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a<b%fc(xyzt)
      END DO
   END FUNCTION real_less_r   

   FUNCTION r_less_eq_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_less_eq_r
      INTEGER xyzt
      r_less_eq_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)<=b%fc(xyzt)
      END DO
   END FUNCTION r_less_eq_r   

   FUNCTION r_less_eq_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_less_eq_real
      INTEGER xyzt
      r_less_eq_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)<=b
      END DO
   END FUNCTION r_less_eq_real   

   FUNCTION real_less_eq_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_less_eq_r
      INTEGER xyzt
      real_less_eq_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a<=b%fc(xyzt)
      END DO
   END FUNCTION real_less_eq_r   

   FUNCTION r_equal_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_equal_r
      INTEGER xyzt
      r_equal_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)==b%fc(xyzt)
      END DO
   END FUNCTION r_equal_r   

   FUNCTION r_equal_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_equal_real
      INTEGER xyzt
      r_equal_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)==b
      END DO
   END FUNCTION r_equal_real   

   FUNCTION real_equal_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_equal_r
      INTEGER xyzt
      real_equal_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a==b%fc(xyzt)
      END DO
   END FUNCTION real_equal_r   

   FUNCTION r_not_equal_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      LOGICAL r_not_equal_r
      INTEGER xyzt
      r_not_equal_r=a%parity==b%parity.AND.a%parity/=-1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)/=b%fc(xyzt)
      END DO
   END FUNCTION r_not_equal_r   

   FUNCTION r_not_equal_real(a,b)
      TYPE(real_field), INTENT(IN) :: a
      REAL(REAL8), INTENT(IN) :: b
      LOGICAL r_not_equal_real
      INTEGER xyzt
      r_not_equal_real=a%parity==0.OR.a%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a%fc(xyzt)/=b
      END DO
   END FUNCTION r_not_equal_real   

   FUNCTION real_not_equal_r(a,b)
      REAL(REAL8), INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      LOGICAL real_not_equal_r
      INTEGER xyzt
      real_not_equal_r=b%parity==0.OR.b%parity==1
      DO xyzt=0,NXYZT2-1
         context(xyzt)=a/=b%fc(xyzt)
      END DO
    END FUNCTION real_not_equal_r   

   FUNCTION r_xor_r(a,b)
      TYPE(real_field), INTENT(IN) :: a,b
      TYPE(real_field) r_xor_r
      INTEGER xyzt
      IF(a%parity==b%parity) THEN
         r_xor_r%parity=a%parity
      ELSE   
         r_xor_r%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
         IF(context(xyzt)) THEN
            r_xor_r%fc(xyzt)=a%fc(xyzt)
         ELSE
            r_xor_r%fc(xyzt)=b%fc(xyzt)
         ENDIF
      END DO
   END FUNCTION r_xor_r

   FUNCTION c_xor_c(a,b)
      TYPE(complex_field), INTENT(IN) :: a,b
      TYPE(complex_field) c_xor_c
      INTEGER xyzt
      IF(a%parity==b%parity) THEN
         c_xor_c%parity=a%parity
      ELSE   
         c_xor_c%parity=-1
      ENDIF
      DO xyzt=0,NXYZT2-1
         IF(context(xyzt)) THEN
            c_xor_c%fc(xyzt)=a%fc(xyzt)
         ELSE
            c_xor_c%fc(xyzt)=b%fc(xyzt)
         ENDIF
      END DO
   END FUNCTION c_xor_c                      

   FUNCTION f_xor_f(a,b)
      TYPE(fermi_field), INTENT(IN) :: a,b
      TYPE(fermi_field) f_xor_f
      INTEGER xyzt,i
      IF(a%parity==b%parity) THEN
         f_xor_f%parity=a%parity
      ELSE   
         f_xor_f%parity=-1
      ENDIF
      DO i=1,4
        DO xyzt=0,NXYZT2-1
          IF(context(xyzt)) THEN
            f_xor_f%fc(1,xyzt,i)=a%fc(1,xyzt,i)
            f_xor_f%fc(2,xyzt,i)=a%fc(2,xyzt,i)
          ELSE
            f_xor_f%fc(1,xyzt,i)=b%fc(1,xyzt,i)
            f_xor_f%fc(2,xyzt,i)=b%fc(2,xyzt,i)
          ENDIF
        END DO
      END DO
   END FUNCTION f_xor_f                        

   FUNCTION g_xor_g(a,b)
      TYPE(gauge_field), INTENT(IN) :: a,b
      TYPE(gauge_field) g_xor_g
      INTEGER xyzt,i,j
      IF(a%parity==b%parity) THEN
         g_xor_g%parity=a%parity
      ELSE   
         g_xor_g%parity=-1
      ENDIF
      IF(a%dir==b%dir) THEN
         g_xor_g%dir=a%dir
      ELSE   
         g_xor_g%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
        IF(context(xyzt)) THEN
          DO j=1,2
          DO i=1,2
              g_xor_g%fc(i,j,xyzt)=a%fc(i,j,xyzt)
          END DO
          END DO
        ELSE
          DO j=1,2
          DO i=1,2
              g_xor_g%fc(i,j,xyzt)=b%fc(i,j,xyzt)
          END DO
          END DO
        ENDIF
      END DO
   END FUNCTION g_xor_g                      

   FUNCTION ge_xor_ge(a,b)
      TYPE(generator_field), INTENT(IN) :: a,b
      TYPE(generator_field) ge_xor_ge
      INTEGER xyzt,i
      IF(a%parity==b%parity) THEN
         ge_xor_ge%parity=a%parity
      ELSE   
         ge_xor_ge%parity=-1
      ENDIF
      IF(a%dir==b%dir) THEN
         ge_xor_ge%dir=a%dir
      ELSE   
         ge_xor_ge%dir=0
      ENDIF
      DO xyzt=0,NXYZT2-1
         IF(context(xyzt)) THEN
            DO i=1,3
              ge_xor_ge%fc(i,xyzt)=a%fc(i,xyzt)
            END DO
         ELSE
            DO i=1,3
              ge_xor_ge%fc(i,xyzt)=b%fc(i,xyzt)
            END DO
         ENDIF
      END DO
   END FUNCTION ge_xor_ge                      

END MODULE




