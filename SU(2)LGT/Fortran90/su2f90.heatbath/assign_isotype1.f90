!  Program qcdf90, module assign_isotype1, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.
 
! su2 version  2006.08.28 by y.k.

MODULE assign_isotype1

USE precisions
USE constants
USE global_module

   IMPLICIT NONE  

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE g_asgn_g
   END INTERFACE

CONTAINS

   SUBROUTINE g_asgn_g(a,b)
      TYPE(gauge_field), INTENT(INOUT) :: a
      TYPE(gauge_field), INTENT(IN) :: b
      INTEGER i,j,xyzt,k,wr,q,p,eo
      COMPLEX(REAL8) aux(2,2)
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=b%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('+')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=a%fc(i,j,xyzt)+b%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('-')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=a%fc(i,j,xyzt)-b%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('*')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         SELECT CASE(assign_spec)
         CASE(0)
            DO xyzt=0,NXYZT2-1
               DO j=1,2
               DO i=1,2
                  aux(i,j)= &
                          a%fc(i,1,xyzt)*b%fc(1,j,xyzt)+ &
                          a%fc(i,2,xyzt)*b%fc(2,j,xyzt)
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
                  aux(i,j)= &
                          b%fc(i,1,xyzt)*a%fc(1,j,xyzt)+ &
                          b%fc(i,2,xyzt)*a%fc(2,j,xyzt)
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

      CASE('/')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         SELECT CASE(assign_spec)
         CASE(0)
            DO xyzt=0,NXYZT2-1
               DO j=1,2
               DO i=1,2
                  aux(i,j)= &
                          a%fc(i,1,xyzt)*CONJG(b%fc(j,1,xyzt))+ &
                          a%fc(i,2,xyzt)*CONJG(b%fc(j,2,xyzt))
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
                  aux(i,j)= &
                          CONJG(b%fc(1,i,xyzt))*a%fc(1,j,xyzt)+ &
                          CONJG(b%fc(2,i,xyzt))*a%fc(2,j,xyzt)
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

      CASE('u')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF
         a%parity=1-b%parity
         IF(a%dir.NE.b%dir) a%dir=0 
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(assign_spec>4.OR.assign_spec<-4.OR.assign_spec==0) THEN
           PRINT '("Illegal value for assign_spec")'
           STOP
         ENDIF
         IF(.NOT.shift_initialized) CALL shift_initialization()
         IF(assign_spec>0) THEN
           IF(b%dir==0) THEN
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,assign_spec)
               q=xyzt
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,1,q))+ &
                  b%fc(k,2,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,2,q))
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)= &
                  u%uc(wr,assign_spec)%fc(j,1,xyzt)*aux(1,i)+ &
                  u%uc(wr,assign_spec)%fc(j,2,xyzt)*aux(2,i)
               END DO
               END DO
             END DO
           ELSE
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,assign_spec)
               q=xyzt_neighbor(xyzt,wr,b%dir)
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,1,q))+ &
                  b%fc(k,2,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,2,q))
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)= &
                  u%uc(wr,assign_spec)%fc(j,1,xyzt)*aux(1,i)+ &
                  u%uc(wr,assign_spec)%fc(j,2,xyzt)*aux(2,i)
               END DO
               END DO
             END DO
           ENDIF
         ELSE
           IF(b%dir==0) THEN
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,4-assign_spec)
               q=p
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*u%uc(wr,-assign_spec)%fc(1,i,q)+ &
                  b%fc(k,2,p)*u%uc(wr,-assign_spec)%fc(2,i,q)
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,j,p))*aux(1,i)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,j,p))*aux(2,i)
               END DO
               END DO
             END DO
           ELSE
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,4-assign_spec)
               q=xyzt_neighbor(p,b%parity,b%dir)
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*u%uc(wr,-assign_spec)%fc(1,i,q)+ &
                  b%fc(k,2,p)*u%uc(wr,-assign_spec)%fc(2,i,q)
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,j,p))*aux(1,i)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,j,p))*aux(2,i)
               END DO
               END DO
             END DO
           ENDIF
         ENDIF

      CASE('U')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF
         a%parity=1-b%parity
         IF(a%dir.NE.b%dir) a%dir=0 
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(assign_spec>4.OR.assign_spec<-4.OR.assign_spec==0) THEN
           PRINT '("Illegal value for assign_spec")'
           STOP
         ENDIF
         IF(.NOT.shift_initialized) CALL shift_initialization()
         IF(assign_spec>0) THEN
           IF(b%dir==0) THEN
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,assign_spec)
               q=xyzt
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,1,q))+ &
                  b%fc(k,2,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,2,q))
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)=a%fc(j,i,xyzt)+ &
                  u%uc(wr,assign_spec)%fc(j,1,xyzt)*aux(1,i)+ &
                  u%uc(wr,assign_spec)%fc(j,2,xyzt)*aux(2,i)
               END DO
               END DO
             END DO
           ELSE
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,assign_spec)
               q=xyzt_neighbor(xyzt,wr,b%dir)
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,1,q))+ &
                  b%fc(k,2,p)*CONJG(u%uc(b%parity,assign_spec)%fc(i,2,q))
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)=a%fc(j,i,xyzt)+ &
                  u%uc(wr,assign_spec)%fc(j,1,xyzt)*aux(1,i)+ &
                  u%uc(wr,assign_spec)%fc(j,2,xyzt)*aux(2,i)
               END DO
               END DO
             END DO
           ENDIF
         ELSE
           IF(b%dir==0) THEN
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,4-assign_spec)
               q=p
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*u%uc(wr,-assign_spec)%fc(1,i,q)+ &
                  b%fc(k,2,p)*u%uc(wr,-assign_spec)%fc(2,i,q)
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)=a%fc(j,i,xyzt)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,j,p))*aux(1,i)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,j,p))*aux(2,i)
               END DO
               END DO
             END DO
           ELSE 
             DO xyzt=0,NXYZT2-1
               p=xyzt_neighbor(xyzt,wr,4-assign_spec)
               q=xyzt_neighbor(p,b%parity,b%dir)
               DO i=1,2
               DO k=1,2
                 aux(k,i)= &
                  b%fc(k,1,p)*u%uc(wr,-assign_spec)%fc(1,i,q)+ &
                  b%fc(k,2,p)*u%uc(wr,-assign_spec)%fc(2,i,q)
               END DO
               END DO
               DO i=1,2
               DO j=1,2
                 a%fc(j,i,xyzt)=a%fc(j,i,xyzt)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,j,p))*aux(1,i)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,j,p))*aux(2,i)
               END DO
               END DO
             END DO
           ENDIF
         ENDIF

      CASE('t')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
            IF(context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=b%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE('f')
         IF(a%parity.NE.b%parity) a%parity=-1
         IF(a%dir.NE.b%dir) a%dir=0
         DO xyzt=0,NXYZT2-1
            IF(.NOT.context(xyzt)) THEN
               DO j=1,2
               DO i=1,2
                  a%fc(i,j,xyzt)=b%fc(i,j,xyzt)
               END DO
               END DO
            ENDIF
         END DO

      CASE('A')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=CONJG(b%fc(j,i,xyzt))
         END DO
         END DO
         END DO

      CASE('I')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)= &
               CMPLX(-AIMAG(b%fc(i,j,xyzt)), &
                     REAL(b%fc(i,j,xyzt),REAL8),REAL8) 
         END DO
         END DO
         END DO

      CASE('M')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=-b%fc(i,j,xyzt)
         END DO
         END DO
         END DO

      CASE('C')
         a%parity=b%parity
         a%dir=b%dir
         DO xyzt=0,NXYZT2-1
         DO j=1,2
         DO i=1,2
            a%fc(i,j,xyzt)=CONJG(b%fc(i,j,xyzt))
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE g_asgn_g

END MODULE

