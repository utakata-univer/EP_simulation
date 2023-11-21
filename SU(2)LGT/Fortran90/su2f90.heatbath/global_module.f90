
!  Program qcdf90, module global_module, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

MODULE global_module

USE precisions

   IMPLICIT NONE

!    INTEGER, PARAMETER :: NX=4,NY=4,NZ=4,NT2=2
    INTEGER, PARAMETER :: NX=8,NY=8,NZ=8,NT2=4
!   INTEGER, PARAMETER :: NX=16,NY=16,NZ=16,NT2=8

   INTEGER, PARAMETER :: NT=2*NT2,NXYZT=NX*NY*NZ*NT,NXYZT2=NX*NY*NZ*NT2
   INTEGER, PARAMETER :: NCGV=9*NXYZT2,NCFV=12*NXYZT2,NRGV=2*NCGV
   INTEGER, PARAMETER :: NRFV=2*NCFV,NRGEV=8*NXYZT2

   TYPE gauge_field
      INTEGER parity
      INTEGER dir
      COMPLEX(REAL8), DIMENSION(2,2,0:NXYZT2-1) :: fc
   END TYPE
 
   TYPE full_gauge_field
      TYPE(gauge_field), DIMENSION(0:1,4) :: uc
   END TYPE

   TYPE fermi_field
      INTEGER parity
      COMPLEX(REAL8), DIMENSION(2,0:NXYZT2-1,4) :: fc
   END TYPE

   TYPE complex_field
      INTEGER parity
      COMPLEX(REAL8), DIMENSION(0:NXYZT2-1) :: fc
   END TYPE

   TYPE real_field
      INTEGER parity
      REAL(REAL8), DIMENSION(0:NXYZT2-1) :: fc
   END TYPE

   TYPE generator_field
      INTEGER parity
      INTEGER dir
      REAL(REAL8), DIMENSION(3,0:NXYZT2-1) :: fc
   END TYPE

   TYPE matrix
      COMPLEX(REAL8), DIMENSION(2,2) :: mc
   END TYPE   
 
! Global variables:

   TYPE(full_gauge_field) u
  
   CHARACTER assign_type
   INTEGER assign_spec
   DATA assign_type,assign_spec/'=',0/

   INTEGER, DIMENSION(0:NX-1,0:NY-1,0:NZ-1,0:NT-1) :: xyzt_index
   INTEGER, DIMENSION(0:NX-1,0:NY-1,0:NZ-1,0:NT-1) :: xyzt_parity
   INTEGER, DIMENSION(0:NXYZT2-1,0:1,4) :: xyzt_cartesian
   INTEGER, DIMENSION(0:NXYZT2-1,0:1,8) :: xyzt_neighbor
   LOGICAL shift_initialized
   DATA shift_initialized/.FALSE./

   LOGICAL, DIMENSION(0:NXYZT2-1) :: context

   INTEGER(LONG) seed_a, seed_b
   INTEGER(LONG), DIMENSION(0:NXYZT2-1) :: seeds


   CONTAINS

   SUBROUTINE shift_initialization()
      INTEGER :: x,y,z,t,xyzt,xyzt_even,xyzt_odd,eo
      INTEGER :: xf,yf,zf,tf,xb,yb,zb,tb
      xyzt_even=0
      xyzt_odd=0
      DO t=0,NT-1
      DO z=0,NZ-1
      DO y=0,NY-1
      DO x=0,NX-1
         IF(IAND(x+y+z+t,1)==0) THEN
            xyzt_parity(x,y,z,t)=0
            xyzt_index(x,y,z,t)=xyzt_even
            xyzt_cartesian(xyzt_even,0,1)=x
            xyzt_cartesian(xyzt_even,0,2)=y
            xyzt_cartesian(xyzt_even,0,3)=z
            xyzt_cartesian(xyzt_even,0,4)=t
            xyzt_even=xyzt_even+1
         ELSE
            xyzt_parity(x,y,z,t)=1
            xyzt_index(x,y,z,t)=xyzt_odd
            xyzt_cartesian(xyzt_odd,1,1)=x
            xyzt_cartesian(xyzt_odd,1,2)=y
            xyzt_cartesian(xyzt_odd,1,3)=z
            xyzt_cartesian(xyzt_odd,1,4)=t
            xyzt_odd=xyzt_odd+1
          ENDIF
      END DO
      END DO
      END DO
      END DO
      DO eo=0,1
      DO xyzt=0,NXYZT2-1
         x=xyzt_cartesian(xyzt,eo,1)
         y=xyzt_cartesian(xyzt,eo,2)
         z=xyzt_cartesian(xyzt,eo,3)
         t=xyzt_cartesian(xyzt,eo,4)
         xf=x+1; IF(xf==NX) xf=0
         yf=y+1; IF(yf==NY) yf=0
         zf=z+1; IF(zf==NZ) zf=0
         tf=t+1; IF(tf==NT) tf=0
         xb=x-1; IF(xb==-1) xb=NX-1
         yb=y-1; IF(yb==-1) yb=NY-1
         zb=z-1; IF(zb==-1) zb=NZ-1
         tb=t-1; IF(tb==-1) tb=NT-1
         xyzt_neighbor(xyzt,eo,1)=xyzt_index(xf,y,z,t)
         xyzt_neighbor(xyzt,eo,2)=xyzt_index(x,yf,z,t)
         xyzt_neighbor(xyzt,eo,3)=xyzt_index(x,y,zf,t)
         xyzt_neighbor(xyzt,eo,4)=xyzt_index(x,y,z,tf)
         xyzt_neighbor(xyzt,eo,5)=xyzt_index(xb,y,z,t)
         xyzt_neighbor(xyzt,eo,6)=xyzt_index(x,yb,z,t)
         xyzt_neighbor(xyzt,eo,7)=xyzt_index(x,y,zb,t)
         xyzt_neighbor(xyzt,eo,8)=xyzt_index(x,y,z,tb)
      END DO
      END DO
      shift_initialized=.TRUE.
   END SUBROUTINE

END MODULE
