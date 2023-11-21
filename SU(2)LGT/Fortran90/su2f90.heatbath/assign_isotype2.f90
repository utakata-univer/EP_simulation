
!  Program qcdf90, module assign_isotype2, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.
 
! su2 version  2006.08.28 by y.k.

MODULE assign_isotype2

USE precisions
USE constants
USE global_module

   IMPLICIT NONE  

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE f_asgn_f
   END INTERFACE

CONTAINS

   SUBROUTINE f_asgn_f(a,b)
      TYPE(fermi_field), INTENT(INOUT) :: a
      TYPE(fermi_field), INTENT(IN) :: b
      COMPLEX(REAL8) aux1(2),aux2(2),acc1(2),acc2(2)
      INTEGER :: xyzt,i,j,k,wr,q,s
      SELECT CASE(assign_type)

      CASE('=')
         a%parity=b%parity
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
             a%fc(i,xyzt,s)=b%fc(i,xyzt,s)
         END DO
         END DO
         END DO

      CASE('+')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
            a%fc(i,xyzt,s)=a%fc(i,xyzt,s)+b%fc(i,xyzt,s)
         END DO
         END DO
         END DO

      CASE('-')
         IF(a%parity.NE.b%parity) a%parity=-1
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
            a%fc(i,xyzt,s)=a%fc(i,xyzt,s)-b%fc(i,xyzt,s)
         END DO
         END DO
         END DO

      CASE('u')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF    
         a%parity=1-b%parity
         IF(a%parity.NE.1-b%parity) a%parity=1-b%parity
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(.NOT.shift_initialized) CALL shift_initialization()
         IF(assign_spec>4.OR.assign_spec<-4.OR.assign_spec==0) THEN
            PRINT '("Illegal value for assign_spec")'
            STOP
         ENDIF
         IF(assign_spec>0) THEN
           DO s=1,4
             DO xyzt=0,NXYZT2-1
               q=xyzt_neighbor(xyzt,wr,assign_spec)
               DO i=1,2
                 a%fc(i,xyzt,s)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*b%fc(1,q,s)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*b%fc(2,q,s)
               END DO
             END DO
           END DO
         ELSE
           DO s=1,4
             DO xyzt=0,NXYZT2-1
               q=xyzt_neighbor(xyzt,wr,4-assign_spec)
               DO i=1,2
                 a%fc(i,xyzt,s)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*b%fc(1,q,s)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*b%fc(2,q,s)
               END DO
             END DO
           END DO
         ENDIF

      CASE('U')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF    
         a%parity=1-b%parity
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(.NOT.shift_initialized) CALL shift_initialization()
         IF(assign_spec>4.OR.assign_spec<-4.OR.assign_spec==0) THEN
            PRINT '("Illegal value for assign_spec")'
            STOP
         ENDIF
         IF(assign_spec>0) THEN
           DO s=1,4
             DO xyzt=0,NXYZT2-1
               q=xyzt_neighbor(xyzt,wr,assign_spec)
               DO i=1,2
                 a%fc(i,xyzt,s)=a%fc(i,xyzt,s)+ &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*b%fc(1,q,s)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*b%fc(2,q,s)
               END DO
             END DO
           END DO
         ELSE
           DO s=1,4
             DO xyzt=0,NXYZT2-1
               q=xyzt_neighbor(xyzt,wr,4-assign_spec)
               DO i=1,2
                 a%fc(i,xyzt,s)=a%fc(i,xyzt,s)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*b%fc(1,q,s)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*b%fc(2,q,s)
               END DO
             END DO
           END DO
         ENDIF

      CASE('w','W')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF
         IF(assign_type=='w') THEN
            DO s=1,4
            DO xyzt=0,NXYZT2-1
            DO i=1,2
              a%fc(i,xyzt,s)=(0._REAL8,0._REAL8)
            END DO
            END DO
            END DO
         ENDIF
         a%parity=1-b%parity
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(.NOT.shift_initialized) CALL shift_initialization()
         SELECT CASE(assign_spec)
         CASE(1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,3)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)-acc2(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)-acc1(i)
            END DO
          END DO
         CASE(2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)-CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+CMPLX(-AIMAG(acc2(i)), &
                                    REAL(acc2(i),REAL8),REAL8)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)-CMPLX(-AIMAG(acc1(i)), &
                                    REAL(acc1(i),REAL8),REAL8)
            END DO
          END DO
         CASE(3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)-acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc2(i)
            END DO
          END DO
         CASE(4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,3)
              aux2(i)=2._REAL8*b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc2(i)
            END DO
          END DO
         CASE(-1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,3)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc2(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc1(i)
            END DO
          END DO
         CASE(-2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)+CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)-CMPLX(-AIMAG(acc2(i)), &
                                    REAL(acc2(i),REAL8),REAL8)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+CMPLX(-AIMAG(acc1(i)), &
                                    REAL(acc1(i),REAL8),REAL8)
            END DO
          END DO
         CASE(-3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)-acc2(i)
            END DO
          END DO
         CASE(-4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,1)
              aux2(i)=2._REAL8*b%fc(i,q,2)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
            END DO
          END DO
         CASE DEFAULT
           PRINT '("Illegal value for assign_spec")'
           STOP
         END SELECT

      CASE('x','X')
         IF(.NOT.(b%parity==0.OR.b%parity==1)) THEN
           PRINT '("Illegal value for parity")'
           STOP
         ENDIF
         IF(assign_type=='x') THEN
            DO s=1,4
            DO xyzt=0,NXYZT2-1
            DO i=1,2
              a%fc(i,xyzt,s)=(0._REAL8,0._REAL8)
            END DO
            END DO
            END DO
         ENDIF
         a%parity=1-b%parity
         wr=1-b%parity  ! wr is set to 1-b%parity 
         IF(.NOT.shift_initialized) CALL shift_initialization()
         SELECT CASE(assign_spec)
         CASE(1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,3)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc2(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc1(i)
            END DO
          END DO
         CASE(2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)+CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3) &
                  -CMPLX(-AIMAG(acc2(i)),REAL(acc2(i),REAL8),REAL8)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4) &
                  +CMPLX(-AIMAG(acc1(i)),REAL(acc1(i),REAL8),REAL8)
            END DO
          END DO
         CASE(3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)-acc2(i)
            END DO
          END DO
         CASE(4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,assign_spec)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,1)
              aux2(i)=2._REAL8*b%fc(i,q,2)
            END DO
            DO i=1,2
              acc1(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux1(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux1(2)
              acc2(i)= &
                      u%uc(wr,assign_spec)%fc(i,1,xyzt)*aux2(1)+ &
                      u%uc(wr,assign_spec)%fc(i,2,xyzt)*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
            END DO
          END DO
         CASE(-1)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,4)
              aux2(i)=b%fc(i,q,2)-b%fc(i,q,3)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)-acc2(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)-acc1(i)
            END DO
          END DO
         CASE(-2)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)+CMPLX(-AIMAG(b%fc(i,q,4)), &
                                    REAL(b%fc(i,q,4),REAL8),REAL8)
              aux2(i)=b%fc(i,q,2)-CMPLX(-AIMAG(b%fc(i,q,3)), &
                                    REAL(b%fc(i,q,3),REAL8),REAL8)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3) &
                  +CMPLX(-AIMAG(acc2(i)),REAL(acc2(i),REAL8),REAL8)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4) &
                  -CMPLX(-AIMAG(acc1(i)),REAL(acc1(i),REAL8),REAL8)
            END DO
          END DO
         CASE(-3)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=b%fc(i,q,1)-b%fc(i,q,3)
              aux2(i)=b%fc(i,q,2)+b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,1)=a%fc(i,xyzt,1)+acc1(i)
              a%fc(i,xyzt,2)=a%fc(i,xyzt,2)+acc2(i)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)-acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc2(i)
            END DO
          END DO
         CASE(-4)
          DO xyzt=0,NXYZT2-1
            q=xyzt_neighbor(xyzt,wr,4-assign_spec)
            DO i=1,2
              aux1(i)=2._REAL8*b%fc(i,q,3)
              aux2(i)=2._REAL8*b%fc(i,q,4)
            END DO
            DO i=1,2
              acc1(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux1(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux1(2)
              acc2(i)= &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(1,i,q))*aux2(1)+ &
                  CONJG(u%uc(b%parity,-assign_spec)%fc(2,i,q))*aux2(2)
              a%fc(i,xyzt,3)=a%fc(i,xyzt,3)+acc1(i)
              a%fc(i,xyzt,4)=a%fc(i,xyzt,4)+acc2(i)
            END DO
          END DO
         CASE DEFAULT
           PRINT '("Illegal value for assign_spec")'
           STOP
         END SELECT

      CASE('C')
         a%parity=b%parity
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
            a%fc(i,xyzt,s)=CONJG(b%fc(i,xyzt,s))
         END DO
         END DO
         END DO

      CASE('I')
         a%parity=b%parity
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)= &
             CMPLX(-AIMAG(b%fc(i,xyzt,s)), &
                     REAL(b%fc(i,xyzt,s),REAL8),REAL8) 
         END DO
         END DO
         END DO

      CASE('M')
         a%parity=b%parity
         DO s=1,4
         DO xyzt=0,NXYZT2-1
         DO i=1,2
           a%fc(i,xyzt,s)=-b%fc(i,xyzt,s)
         END DO
         END DO
         END DO

      CASE DEFAULT
         PRINT '("Illegal value for assign_type")'
         STOP
      END SELECT
      assign_type='='
      assign_spec=0
   END SUBROUTINE f_asgn_f

END MODULE

