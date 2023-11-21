
!  Program qcdf90, module constants, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

MODULE constants

USE precisions

   IMPLICIT NONE  

   REAL(REAL8), PARAMETER :: PI=3.141592653589793_REAL8,               &
      TWOPI=6.283185307179586_REAL8, PI2=1.570796326794896_REAL8,      &
      SQRT2=1.414213562373095_REAL8, SQRT22=0.7071067811865475_REAL8,  &
      SQRT3=1.732050807568877_REAL8, SQRT33=0.5773502691896257_REAL8,  &
      TWOSQRT33=1.154700538379251_REAL8

   COMPLEX(REAL8), PARAMETER :: IU=(0._REAL8,1._REAL8)

   COMPLEX(REAL8), DIMENSION(2,2), PARAMETER ::                      &

      ZERO_m=RESHAPE(SOURCE=(/&

      (0._REAL8,0._REAL8),(0._REAL8,0._REAL8),   &
      (0._REAL8,0._REAL8),(0._REAL8,0._REAL8)/), &

      SHAPE=(/2,2/)),                                                &

      UNIT=RESHAPE(SOURCE=(/&

      (1._REAL8,0._REAL8),(0._REAL8,0._REAL8),   &
      (0._REAL8,0._REAL8),(1._REAL8,0._REAL8)/), &

      SHAPE=(/2,2/))

   COMPLEX(REAL8), DIMENSION(2), PARAMETER :: ZERO_v=(/&        
      (0._REAL8,0._REAL8),(0._REAL8,0._REAL8)/)

   REAL(REAL8), DIMENSION(3), PARAMETER :: ZERO_ge=(/&
       0._REAL8,0._REAL8,0._REAL8/)

   COMPLEX(REAL8), DIMENSION(2,2,3), PARAMETER ::                       &

      LAMBDA=RESHAPE(SOURCE=(/&

      (0._REAL8, 0._REAL8), (1._REAL8, 0._REAL8), &
      (1._REAL8, 0._REAL8), (0._REAL8, 0._REAL8), &

      (0._REAL8, 0._REAL8), (0._REAL8, 1._REAL8), &
      (0._REAL8,-1._REAL8), (0._REAL8, 0._REAL8), &

      (1._REAL8, 0._REAL8), (0._REAL8, 0._REAL8), &
      (0._REAL8, 0._REAL8), (-1._REAL8, 0._REAL8)/), &

      SHAPE=(/2,2,3/))

   COMPLEX(REAL8), DIMENSION(4,4,5), PARAMETER ::                       &

      GAMMA=RESHAPE(SOURCE=(/&

      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (1._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(1._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&

      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 1._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8,-1._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 1._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8,-1._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&

      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(1._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                  (-1._REAL8, 0._REAL8),&
      (1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(-1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),  &
                                                   (0._REAL8, 0._REAL8),&

      (1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(1._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(-1._REAL8, 0._REAL8),  &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                  (-1._REAL8, 0._REAL8),&

      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 1._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 1._REAL8),&
      (0._REAL8,-1._REAL8),(0._REAL8, 0._REAL8),(0._REAL8, 0._REAL8),   &
                                                   (0._REAL8, 0._REAL8),&
      (0._REAL8, 0._REAL8),(0._REAL8,-1._REAL8),(0._REAL8, 0._REAL8),   &
                                                (0._REAL8, 0._REAL8)/), &

      SHAPE=(/4,4,5/))

   COMPLEX(REAL8), DIMENSION(2,2), PARAMETER ::                      &

      SIGMA1=RESHAPE(SOURCE=(/&

      (0._REAL8,0._REAL8),(1._REAL8,0._REAL8),   &
      (1._REAL8,0._REAL8),(0._REAL8,0._REAL8)/), &

      SHAPE=(/2,2/)),                                                &

      SIGMA2=RESHAPE(SOURCE=(/&

      (0._REAL8,0._REAL8),(0._REAL8,-1._REAL8),   &
      (0._REAL8,1._REAL8),(0._REAL8,0._REAL8)/), &

      SHAPE=(/2,2/)),                                                &

      SIGMA3=RESHAPE(SOURCE=(/&

      (1._REAL8,0._REAL8),(0._REAL8,0._REAL8),   &
      (0._REAL8,0._REAL8),(-1._REAL8,0._REAL8)/), &

      SHAPE=(/2,2/))

END MODULE








