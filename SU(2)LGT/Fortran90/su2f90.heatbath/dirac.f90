!  Program qcdf90, module dirac, version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

 MODULE dirac

   USE precisions
   USE global_module
   USE field_algebra
   USE assign_isotype2

   IMPLICIT NONE

   INTERFACE OPERATOR (.Dirac.)
      MODULE PROCEDURE dirac_f
   END INTERFACE

   INTERFACE OPERATOR (.Xdirac.)
      MODULE PROCEDURE g5_dirac_g5_f
   END INTERFACE


   CONTAINS

   FUNCTION dirac_f(a)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(fermi_field) :: dirac_f
      INTEGER xyzt,m
      assign_type='w';assign_spec=-4; dirac_f=a
      DO m=-3,4
        IF(m.NE.0) THEN
          assign_type='W';assign_spec=m; dirac_f=a
        ENDIF
      END DO
      dirac_f%parity=1-a%parity
   END FUNCTION dirac_f

   FUNCTION g5_dirac_g5_f(a)
      TYPE(fermi_field), INTENT(IN) :: a
      TYPE(fermi_field) :: g5_dirac_g5_f
      INTEGER xyzt,m
      assign_type='x';assign_spec=-4; g5_dirac_g5_f=a
      DO m=-3,4
        IF(m.NE.0) THEN
          assign_type='X';assign_spec=m; g5_dirac_g5_f=a
        ENDIF
      END DO
      g5_dirac_g5_f%parity=1-a%parity
   END FUNCTION g5_dirac_g5_f

END MODULE

