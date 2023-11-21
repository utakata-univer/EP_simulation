
!  Program qcdf90, module random_numbers, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.
 
MODULE random_numbers

USE precisions
USE global_module
USE constants

   IMPLICIT NONE

                    ! Unary operator seed: when acting on seed (long integer)
                    ! sets the global variable seeds to the rand48
                    ! sequence beginning with seed and also sets the
                    ! global variables seed_a, seed_b to the appropriate
                    ! multiplier and constant term needed to produce
                    ! skips by NX*NY*NZ*NT/2 in the sequence of pseudorandom
                    ! numbers.  It returns .TRUE.
                    ! When acting on a logical variable set to .TRUE. it
                    ! returns the current seed (=seeds(0,0,0,0)) that must
                    ! be used to restart the sequence of random numbers.
                    ! If the argument is .FALSE. it returns 0.  

   INTERFACE OPERATOR(.Seed.)
      MODULE PROCEDURE seed_in, seed_out
   END INTERFACE

                    ! Unary operator .Rand.: returns a real field of
                    ! pseudorandom numbers with uniform distribution
                    ! between 0 and the argument passed to .Rand. and
                    ! upgrades the global variable seeds.  
                    ! Unary operator .Gauss.: returns a real field of
                    ! pseudorandom numbers with gaussian distribution
                    ! of mean square deviation equal to the argument of
                    ! .Gauss. and upgrades the global variable seeds.  
                    ! The argument can be a real constant, in which case
                    ! the parity of the result is undefined, or a real
                    ! field, in which case the parity of the result
                    ! is the same as the parity of the argument.
                    ! Note that .Seed. must be used to initialize
                    ! the generation of psuedorandom numbers prior
                    ! to the use of .Rand. or .Gauss.
                    ! .Ggauss. works like .Gauss. but fills with
                    ! gaussian random numbers the components of a
                    ! generator_field, setting its direction equal to 0.


   INTERFACE OPERATOR(.Rand.)
      MODULE PROCEDURE uniform_real, uniform_r
   END INTERFACE
   INTERFACE OPERATOR(.Gauss.)
      MODULE PROCEDURE gaussian_real, gaussian_r
   END INTERFACE
   INTERFACE OPERATOR(.Ggauss.)
      MODULE PROCEDURE ggaussian_real, ggaussian_r
   END INTERFACE

CONTAINS

   FUNCTION seed_in(saved_seed)
      INTEGER(LONG), INTENT(IN) :: saved_seed
      LOGICAL seed_in
      INTEGER(LONG) :: a,b,c
      INTEGER aux,i,m,xyzt

      seed_in=.FALSE.

      a=25214903917_LONG
      b=11_LONG
      c=saved_seed

      DO xyzt=0,NXYZT2-1
         seeds(xyzt)=c
         c=IAND(a*c+b,281474976710655_LONG)
      END DO

      seed_a=1
      seed_b=0
      aux=NXYZT2

      DO i=1,48
         m=IAND(aux,1)
         aux=ISHFT(aux,-1)
         IF(m.EQ.1) THEN
            seed_b=IAND(a*seed_b+b,281474976710655_LONG)
            seed_a=IAND(a*seed_a,281474976710655_LONG)
         ENDIF    
         IF(aux.EQ.0) THEN
            GOTO 10
         ENDIF
         b=IAND(a*b+b,281474976710655_LONG)
         a=IAND(a*a,281474976710655_LONG)
      ENDDO

 10   seed_in=.TRUE.

   END FUNCTION seed_in
 
   FUNCTION seed_out(a)
      LOGICAL, INTENT(IN) :: a
      INTEGER(LONG) seed_out
         IF(a) THEN
         seed_out=seeds(0)
      ELSE
         seed_out=0
      ENDIF
   END FUNCTION seed_out

   FUNCTION uniform_real(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      REAL(REAL8), INTENT(IN) :: r
      REAL(REAL8) sf
      TYPE(real_field) uniform_real
      INTEGER xyzt
      sf=r*TWONEG48
      uniform_real%parity=-1
      DO xyzt=0,NXYZT2-1
         seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                 seed_b,281474976710655_LONG)
         uniform_real%fc(xyzt)=sf*seeds(xyzt)
      END DO
   END FUNCTION uniform_real

   FUNCTION uniform_r(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      TYPE(real_field), INTENT(IN) :: r
      TYPE(real_field) uniform_r
      INTEGER xyzt
      uniform_r%parity=r%parity
      DO xyzt=0,NXYZT2-1
         seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                 seed_b,281474976710655_LONG)
         uniform_r%fc(xyzt)=TWONEG48*r%fc(xyzt)*seeds(xyzt)
      END DO
   END FUNCTION uniform_r

   FUNCTION gaussian_real(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      REAL(REAL8), INTENT(IN) :: r
      REAL(REAL8) a,c,s,sf 
      TYPE(real_field) :: gaussian_real
      INTEGER xyzt
      sf=TWOPI*TWONEG48
      gaussian_real%parity=-1
      DO xyzt=0,NXYZT2-1,2
         seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                 seed_b,281474976710655_LONG)
         seeds(xyzt+1)=IAND(seed_a*seeds(xyzt+1)+          &
                                 seed_b,281474976710655_LONG)
         a=r*SQRT(-2._REAL8*LOG(TWONEG48*seeds(xyzt)))
         c=COS(sf*seeds(xyzt+1))
         s=SIN(sf*seeds(xyzt+1))
         gaussian_real%fc(xyzt)=a*c
         gaussian_real%fc(xyzt+1)=a*s
      END DO
   END FUNCTION gaussian_real

   FUNCTION gaussian_r(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      TYPE(real_field), INTENT(IN) :: r
      REAL(REAL8) a,c,s,sf
      TYPE(real_field) :: gaussian_r
      INTEGER xyzt
      sf=TWOPI*TWONEG48
      gaussian_r%parity=r%parity
      DO xyzt=0,NXYZT2-1,2
         seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                 seed_b,281474976710655_LONG)
         seeds(xyzt+1)=IAND(seed_a*seeds(xyzt+1)+          &
                                 seed_b,281474976710655_LONG)
         a=SQRT(-2._REAL8*LOG(TWONEG48*seeds(xyzt)))
         c=COS(sf*seeds(xyzt+1))
         s=SIN(sf*seeds(xyzt+1))
         gaussian_r%fc(xyzt)=a*c*r%fc(xyzt)
         gaussian_r%fc(xyzt+1)=a*s*r%fc(xyzt+1)
      END DO
   END FUNCTION gaussian_r

   FUNCTION ggaussian_real(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      REAL(REAL8), INTENT(IN) :: r
      REAL(REAL8) a,c,s,sf 
      TYPE(generator_field) :: ggaussian_real
      INTEGER xyzt,i
      sf=TWOPI*TWONEG48
      ggaussian_real%parity=-1
      ggaussian_real%dir=0
      DO i=1,3
         DO xyzt=0,NXYZT2-1,2
            seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                    seed_b,281474976710655_LONG)
            seeds(xyzt+1)=IAND(seed_a*seeds(xyzt+1)+          &
                                    seed_b,281474976710655_LONG)
            a=r*SQRT(-2._REAL8*LOG(TWONEG48*seeds(xyzt)))
            c=COS(sf*seeds(xyzt+1))
            s=SIN(sf*seeds(xyzt+1))
            ggaussian_real%fc(i,xyzt)=a*c
            ggaussian_real%fc(i,xyzt+1)=a*s
         END DO
      END DO   
   
   END FUNCTION ggaussian_real

   FUNCTION ggaussian_r(r)
      REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
      TYPE(real_field), INTENT(IN) :: r
      REAL(REAL8) a,c,s,sf 
      TYPE(generator_field) :: ggaussian_r
      INTEGER xyzt,i
      sf=TWOPI*TWONEG48
      ggaussian_r%parity=r%parity
      ggaussian_r%dir=0
      DO i=1,3
         DO xyzt=0,NXYZT2-1,2
            seeds(xyzt)=IAND(seed_a*seeds(xyzt)+          &
                                    seed_b,281474976710655_LONG)
            seeds(xyzt+1)=IAND(seed_a*seeds(xyzt+1)+          &
                                    seed_b,281474976710655_LONG)
            a=SQRT(-2._REAL8*LOG(TWONEG48*seeds(xyzt)))
            c=COS(sf*seeds(xyzt+1))
            s=SIN(sf*seeds(xyzt+1))
            ggaussian_r%fc(i,xyzt)=a*c*r%fc(xyzt)
            ggaussian_r%fc(i,xyzt+1)=a*s*r%fc(xyzt+1)
         END DO
      END DO   
   END FUNCTION ggaussian_r

END MODULE



