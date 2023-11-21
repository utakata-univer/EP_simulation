
!  Program qcdf90, module shift, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.
 
 MODULE shift

   USE precisions
   USE global_module
   IMPLICIT NONE

   INTERFACE OPERATOR (.Cshift.)
      MODULE PROCEDURE cshift_g, cshift_f, cshift_c, &
                       cshift_r, cshift_ge
   END INTERFACE

   INTERFACE OPERATOR (.Ushift.)
      MODULE PROCEDURE ushift_g, ushift_f, cshift_c, &
                       cshift_r
   END INTERFACE

   INTERFACE OPERATOR (.Wshift.)
      MODULE PROCEDURE wshift_f
   END INTERFACE

   INTERFACE OPERATOR (.Xshift.)
      MODULE PROCEDURE g5_wshift_g5_f
   END INTERFACE

   CONTAINS

   FUNCTION Cshift_g(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) :: cshift_g
      INTEGER :: xyzt,q,mu,i,j
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Cshift parameters")'
        STOP
      ENDIF
      cshift_g%parity=1-b%parity
      cshift_g%dir=b%dir
      mu=a; IF(a<0) mu=4-a
      DO xyzt=0,NXYZT2-1
         q=xyzt_neighbor(xyzt,cshift_g%parity,mu)
         DO j=1,2
         DO i=1,2
            cshift_g%fc(i,j,xyzt)=b%fc(i,j,q)
         END DO
         END DO
      END DO
   END FUNCTION cshift_g

   FUNCTION cshift_f(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) :: cshift_f
      INTEGER :: xyzt,q,mu,i,s
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Cshift parameters")'
        STOP
      ENDIF
      cshift_f%parity=1-b%parity
      mu=a; IF(a<0) mu=4-a
      DO s=1,4
         DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,cshift_f%parity,mu)
            DO i=1,2
               cshift_f%fc(i,xyzt,s)=b%fc(i,q,s)
            END DO
         END DO
      END DO
   END FUNCTION cshift_f

   FUNCTION cshift_c(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(complex_field), INTENT(IN) :: b
      TYPE(complex_field) :: cshift_c
      INTEGER :: xyzt,q,mu
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Cshift parameters")'
        STOP
      ENDIF
      cshift_c%parity=1-b%parity
      mu=a; IF(a<0) mu=4-a
      DO xyzt=0,NXYZT2-1
         q=xyzt_neighbor(xyzt,cshift_c%parity,mu)
         cshift_c%fc(xyzt)=b%fc(q)
      END DO
   END FUNCTION cshift_c

   FUNCTION cshift_r(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(real_field), INTENT(IN) :: b
      TYPE(real_field) :: cshift_r
      INTEGER :: xyzt,q,mu
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Cshift parameters")'
        STOP
      ENDIF
      cshift_r%parity=1-b%parity
      mu=a; IF(a<0) mu=4-a
      DO xyzt=0,NXYZT2-1
         q=xyzt_neighbor(xyzt,cshift_r%parity,mu)
         cshift_r%fc(xyzt)=b%fc(q)
      END DO
   END FUNCTION cshift_r

   FUNCTION cshift_ge(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(generator_field), INTENT(IN) :: b
      TYPE(generator_field) :: cshift_ge
      INTEGER :: xyzt,q,mu,i
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Cshift parameters")'
        STOP
      ENDIF
      cshift_ge%parity=1-b%parity
      cshift_ge%dir=b%dir
      mu=a; IF(a<0) mu=4-a
      DO xyzt=0,NXYZT2-1
         q=xyzt_neighbor(xyzt,cshift_ge%parity,mu)
         DO i=1,3
            cshift_ge%fc(i,xyzt)=b%fc(i,q)
         END DO
      END DO
   END FUNCTION cshift_ge

   FUNCTION ushift_g(a,b)
      INTEGER, INTENT(IN) :: a  
      TYPE(gauge_field), INTENT(IN) :: b
      TYPE(gauge_field) ushift_g
      COMPLEX(REAL8) temp(2,2)
      INTEGER :: xyzt,q,p,i,j,k,wr
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Ushift parameters")'
        STOP
      ENDIF
      ushift_g%dir=b%dir
      ushift_g%parity=1-b%parity
      wr=ushift_g%parity  ! wr is set to 1-b%parity 
      IF(a>0) THEN
        IF(b%dir==0) THEN
          DO xyzt=0,NXYZT2-1
            p=xyzt_neighbor(xyzt,wr,a)
            q=xyzt
            DO i=1,2
            DO k=1,2
              temp(k,i)= &
               b%fc(k,1,p)*CONJG(u%uc(b%parity,a)%fc(i,1,q))+ &
               b%fc(k,2,p)*CONJG(u%uc(b%parity,a)%fc(i,2,q))
            END DO
            END DO
            DO i=1,2
            DO j=1,2
              ushift_g%fc(j,i,xyzt)= &
                u%uc(wr,a)%fc(j,1,xyzt)*temp(1,i)+ &
                u%uc(wr,a)%fc(j,2,xyzt)*temp(2,i)
            END DO
            END DO
          END DO
        ELSE
          DO xyzt=0,NXYZT2-1
            p=xyzt_neighbor(xyzt,wr,a)
            q=xyzt_neighbor(xyzt,wr,b%dir)
            DO i=1,2
            DO k=1,2
              temp(k,i)= &
               b%fc(k,1,p)*CONJG(u%uc(b%parity,a)%fc(i,1,q))+ &
               b%fc(k,2,p)*CONJG(u%uc(b%parity,a)%fc(i,2,q))
            END DO
            END DO
            DO i=1,2
            DO j=1,2
              ushift_g%fc(j,i,xyzt)= &
                u%uc(wr,a)%fc(j,1,xyzt)*temp(1,i)+ &
                u%uc(wr,a)%fc(j,2,xyzt)*temp(2,i)
            END DO
            END DO
          END DO
        ENDIF
      ELSE
        IF(b%dir==0) THEN
          DO xyzt=0,NXYZT2-1
            p=xyzt_neighbor(xyzt,wr,4-a)
            q=p
            DO i=1,2
            DO k=1,2
              temp(k,i)= &
               b%fc(k,1,p)*u%uc(wr,-a)%fc(1,i,q)+ &
               b%fc(k,2,p)*u%uc(wr,-a)%fc(2,i,q)
            END DO
            END DO
            DO i=1,2
            DO j=1,2
              ushift_g%fc(j,i,xyzt)= &
                CONJG(u%uc(b%parity,-a)%fc(1,j,p))*temp(1,i)+ &
                CONJG(u%uc(b%parity,-a)%fc(2,j,p))*temp(2,i)
            END DO
            END DO
          END DO
        ELSE
          DO xyzt=0,NXYZT2-1
            p=xyzt_neighbor(xyzt,wr,4-a)
            q=xyzt_neighbor(p,b%parity,b%dir)
            DO i=1,2
            DO k=1,2
              temp(k,i)= &
               b%fc(k,1,p)*u%uc(wr,-a)%fc(1,i,q)+ &
               b%fc(k,2,p)*u%uc(wr,-a)%fc(2,i,q)
            END DO
            END DO
            DO i=1,2
            DO j=1,2
              ushift_g%fc(j,i,xyzt)= &
                CONJG(u%uc(b%parity,-a)%fc(1,j,p))*temp(1,i)+ &
                CONJG(u%uc(b%parity,-a)%fc(2,j,p))*temp(2,i)
            END DO
            END DO
          END DO
        ENDIF
      ENDIF
   END FUNCTION ushift_g

   FUNCTION ushift_f(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) :: ushift_f
      INTEGER :: xyzt,q,i,s,wr
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Ushift parameters")'
        STOP
      ENDIF
      ushift_f%parity=1-b%parity
      wr=ushift_f%parity  !wr is set to 1-b%parity 
      IF(a>0) THEN
        DO s=1,4
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              ushift_f%fc(i,xyzt,s)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*b%fc(1,q,s)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*b%fc(2,q,s)
            END DO
          END DO
        END DO
      ELSE
        DO s=1,4
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              ushift_f%fc(i,xyzt,s)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*b%fc(1,q,s)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*b%fc(2,q,s)
            END DO
          END DO
        END DO
      ENDIF
   END FUNCTION ushift_f

   FUNCTION wshift_f(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) :: wshift_f
      COMPLEX(REAL8), DIMENSION(2) :: aux1,aux2
      INTEGER xyzt,q,i,wr
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Wshift parameters")'
        STOP
      ENDIF
      wshift_f%parity=1-b%parity
      wr=wshift_f%parity  !wr is set to 1-b%parity 
      SELECT CASE (a)
      CASE(1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,3)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              wshift_f%fc(i,xyzt,3)=-wshift_f%fc(i,xyzt,2)
              wshift_f%fc(i,xyzt,4)=-wshift_f%fc(i,xyzt,1)
            END DO
          END DO
      CASE(2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)-CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              wshift_f%fc(i,xyzt,3)=CMPLX(-AIMAG(wshift_f%fc(i,xyzt,2)), &
                                    REAL(wshift_f%fc(i,xyzt,2),REAL8),REAL8)
              wshift_f%fc(i,xyzt,4)=-CMPLX(-AIMAG(wshift_f%fc(i,xyzt,1)), &
                                    REAL(wshift_f%fc(i,xyzt,1),REAL8),REAL8)
            END DO
          END DO
      CASE(3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,4)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              wshift_f%fc(i,xyzt,3)=-wshift_f%fc(i,xyzt,1)
              wshift_f%fc(i,xyzt,4)=wshift_f%fc(i,xyzt,2)
            END DO
          END DO
      CASE(4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,3)
              aux2(i)=2._REAL8*b%fc(i,q,4)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,3)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              wshift_f%fc(i,xyzt,4)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              wshift_f%fc(i,xyzt,1)=(0._REAL8,0._REAL8)
              wshift_f%fc(i,xyzt,2)=(0._REAL8,0._REAL8)
            END DO
          END DO
      CASE(-1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,3)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              wshift_f%fc(i,xyzt,3)=wshift_f%fc(i,xyzt,2)
              wshift_f%fc(i,xyzt,4)=wshift_f%fc(i,xyzt,1)
            END DO
          END DO
      CASE(-2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)+CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              wshift_f%fc(i,xyzt,3)=-CMPLX(-AIMAG(wshift_f%fc(i,xyzt,2)), &
                                    REAL(wshift_f%fc(i,xyzt,2),REAL8),REAL8)
              wshift_f%fc(i,xyzt,4)=CMPLX(-AIMAG(wshift_f%fc(i,xyzt,1)), &
                                    REAL(wshift_f%fc(i,xyzt,1),REAL8),REAL8)
            END DO
          END DO
      CASE(-3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,4)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              wshift_f%fc(i,xyzt,3)=wshift_f%fc(i,xyzt,1)
              wshift_f%fc(i,xyzt,4)=-wshift_f%fc(i,xyzt,2)
            END DO
          END DO
      CASE(-4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,1)
              aux2(i)=2._REAL8*b%fc(i,q,2)
            END DO
            DO i=1,2
              wshift_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              wshift_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              wshift_f%fc(i,xyzt,3)=(0._REAL8,0._REAL8)
              wshift_f%fc(i,xyzt,4)=(0._REAL8,0._REAL8)
            END DO
          END DO
       CASE DEFAULT
           PRINT *,'shift variable a out of range'
           STOP
       END SELECT    
   END FUNCTION wshift_f

   FUNCTION g5_wshift_g5_f(a,b)
      INTEGER, INTENT(IN) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      TYPE(fermi_field) :: g5_wshift_g5_f
      COMPLEX(REAL8), DIMENSION(2) :: aux1,aux2
      INTEGER xyzt,q,i,s,wr
      IF(.NOT.shift_initialized) CALL shift_initialization()
      IF((.NOT.(b%parity==0.OR.b%parity==1)).OR.a<-4.OR.a>4.OR.a==0) THEN
        PRINT '("Wrong Xshift parameters")'
        STOP
      ENDIF
      g5_wshift_g5_f%parity=1-b%parity
      wr=g5_wshift_g5_f%parity  !wr is set to 1-b%parity 
      SELECT CASE (a)
      CASE(1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,3)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)=g5_wshift_g5_f%fc(i,xyzt,2)
              g5_wshift_g5_f%fc(i,xyzt,4)=g5_wshift_g5_f%fc(i,xyzt,1)
            END DO
          END DO
      CASE(2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)+CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)= &
                  -CMPLX(-AIMAG(g5_wshift_g5_f%fc(i,xyzt,2)), &
                  REAL(g5_wshift_g5_f%fc(i,xyzt,2),REAL8),REAL8)
              g5_wshift_g5_f%fc(i,xyzt,4)= &
                  CMPLX(-AIMAG(g5_wshift_g5_f%fc(i,xyzt,1)), &
                  REAL(g5_wshift_g5_f%fc(i,xyzt,1),REAL8),REAL8)
            END DO
          END DO
      CASE(3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,4)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)=g5_wshift_g5_f%fc(i,xyzt,1)
              g5_wshift_g5_f%fc(i,xyzt,4)=-g5_wshift_g5_f%fc(i,xyzt,2)
            END DO
          END DO
      CASE(4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,a)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,1)
              aux2(i)=2._REAL8*b%fc(i,q,2)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                      u%uc(wr,a)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,a)%fc(i,2,xyzt)*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)=(0._REAL8,0._REAL8)
              g5_wshift_g5_f%fc(i,xyzt,4)=(0._REAL8,0._REAL8)
            END DO
          END DO
      CASE(-1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,3)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)=-g5_wshift_g5_f%fc(i,xyzt,2)
              g5_wshift_g5_f%fc(i,xyzt,4)=-g5_wshift_g5_f%fc(i,xyzt,1)
            END DO
          END DO
      CASE(-2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)-CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)= &
                  CMPLX(-AIMAG(g5_wshift_g5_f%fc(i,xyzt,2)), &
                  REAL(g5_wshift_g5_f%fc(i,xyzt,2),REAL8),REAL8)
              g5_wshift_g5_f%fc(i,xyzt,4)= &
                  -CMPLX(-AIMAG(g5_wshift_g5_f%fc(i,xyzt,1)), &
                  REAL(g5_wshift_g5_f%fc(i,xyzt,1),REAL8),REAL8)
            END DO
          END DO
      CASE(-3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,4)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,1)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,2)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,3)=-g5_wshift_g5_f%fc(i,xyzt,1)
              g5_wshift_g5_f%fc(i,xyzt,4)=g5_wshift_g5_f%fc(i,xyzt,2)
            END DO
          END DO
      CASE(-4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-a)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,3)
              aux2(i)=2._REAL8*b%fc(i,q,4)
            END DO
            DO i=1,2
              g5_wshift_g5_f%fc(i,xyzt,3)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux1(2)
              g5_wshift_g5_f%fc(i,xyzt,4)= &
                  CONJG(u%uc(b%parity,-a)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-a)%fc(2,i,q))*aux2(2)
              g5_wshift_g5_f%fc(i,xyzt,1)=(0._REAL8,0._REAL8)
              g5_wshift_g5_f%fc(i,xyzt,2)=(0._REAL8,0._REAL8)
            END DO
          END DO
       CASE DEFAULT
           PRINT *,'shift variable a out of range'
           STOP
       END SELECT    
  END FUNCTION g5_wshift_g5_f

END MODULE

