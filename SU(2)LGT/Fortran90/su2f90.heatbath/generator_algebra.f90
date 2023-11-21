
!  Program qcdf90, module generator_algebra, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

MODULE generator_algebra

USE precisions
USE constants
USE global_module

   IMPLICIT NONE  


   INTERFACE OPERATOR(.Matrix.)
      MODULE PROCEDURE g_from_ge
   END INTERFACE

   INTERFACE OPERATOR(.Generator.)
      MODULE PROCEDURE ge_from_g
   END INTERFACE

   INTERFACE OPERATOR(.Sq.)
      MODULE PROCEDURE sq_ge
   END INTERFACE

   INTERFACE OPERATOR(.Exp.)
      MODULE PROCEDURE exp_ge
   END INTERFACE

CONTAINS

   FUNCTION g_from_ge(a)
      TYPE(generator_field), INTENT(IN) :: a
      TYPE(gauge_field) g_from_ge
      INTEGER xyzt
      g_from_ge%parity=a%parity    
      g_from_ge%dir=a%dir            
      DO xyzt=0,NXYZT2-1
        g_from_ge%fc(1,1,xyzt)=a%fc(3,xyzt)
        g_from_ge%fc(1,2,xyzt)=CMPLX(a%fc(1,xyzt),-a%fc(2,xyzt),REAL8)
        g_from_ge%fc(2,1,xyzt)=CMPLX(a%fc(1,xyzt), a%fc(2,xyzt),REAL8)
        g_from_ge%fc(2,2,xyzt)=-a%fc(3,xyzt)
      END DO
   END FUNCTION g_from_ge

   FUNCTION ge_from_g(a)
      TYPE(gauge_field), INTENT(IN) :: a
      TYPE(generator_field) ge_from_g
      INTEGER xyzt
      ge_from_g%parity=a%parity    
      ge_from_g%dir=a%dir
      DO xyzt=0,NXYZT2-1
         ge_from_g%fc(1,xyzt)=0.5_REAL8*AIMAG(a%fc(1,2,xyzt)  &
            +a%fc(2,1,xyzt))
         ge_from_g%fc(2,xyzt)=0.5_REAL8*REAL(a%fc(1,2,xyzt)   &
            -a%fc(2,1,xyzt),REAL8)
         ge_from_g%fc(3,xyzt)=0.5_REAL8*AIMAG(a%fc(1,1,xyzt)  &
            -a%fc(2,2,xyzt))
      END DO
   END FUNCTION ge_from_g

   FUNCTION sq_ge(a)
      TYPE(generator_field), INTENT(IN) :: a
      TYPE(generator_field) sq_ge
      INTEGER xyzt
      sq_ge%parity=a%parity    
      sq_ge%dir=a%dir            
      DO xyzt=0,NXYZT2-1
         sq_ge%fc(1,xyzt)=0.0_REAL8
         sq_ge%fc(2,xyzt)=0.0_REAL8
         sq_ge%fc(3,xyzt)=0.0_REAL8
      END DO         
   END FUNCTION sq_ge

   FUNCTION exp_ge(h)
      TYPE(generator_field), INTENT(IN) :: h
      TYPE(gauge_field) exp_ge
      REAL(REAL8) s,c1,c2
      INTEGER xyzt,i,j,k
      COMPLEX(REAL8), DIMENSION(2,2) :: ms

      exp_ge%parity=h%parity    
      exp_ge%dir=h%dir            

      DO xyzt=0,NXYZT2-1

         s=h%fc(1,xyzt)**2
         DO i=2,3
            s=s+h%fc(i,xyzt)**2
         END DO
         s=SQRT(s)
     
         c1=COS(s)
         IF(ABS(s)>0.00000001_REAL8) THEN
            c2=SIN(s)/s
         ELSE
            c2=1._REAL8
         ENDIF 
     
!  ms=c1*UNIT+IU*c2*(.Matrix.hs.), inlined:
     

         exp_ge%fc(1,1,xyzt)=c1+c2*CMPLX(0._REAL8,  h%fc(3,xyzt),REAL8)
         exp_ge%fc(2,2,xyzt)=c1+c2*CMPLX(0._REAL8, -h%fc(3,xyzt),REAL8)
         exp_ge%fc(1,2,xyzt)=c2*CMPLX( h%fc(2,xyzt),h%fc(1,xyzt),REAL8)
         exp_ge%fc(2,1,xyzt)=c2*CMPLX(-h%fc(2,xyzt),h%fc(1,xyzt),REAL8)

      END DO         
   END FUNCTION exp_ge

END MODULE

