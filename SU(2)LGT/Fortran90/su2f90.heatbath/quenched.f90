!  Program Qcdf90_quenched version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  May 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.

PROGRAM Qcdf90_quenched

   USE precisions
   USE constants
   USE global_module
   USE field_algebra
   USE generator_algebra
   USE conditionals
   USE shift
   USE random_numbers
   USE assign_mixed
   USE assign_isotype1

   IMPLICIT NONE
   TYPE(gauge_field):: staple,g_old,g_new
   TYPE(real_field):: plaq_old,plaq_new,bf_ratio,rand
   TYPE(generator_field):: ge
   TYPE(matrix) :: zero_matrix, unit_matrix
   LOGICAL l_test,l_seed
   REAL(REAL8) clock_dcl,clock_upd,clock_plaq
   REAL(REAL8) beta,saved_beta,hp,av_plaq,aux,range_small,range_unit
   CHARACTER(LEN=64) in_filename,out_filename
   CHARACTER(LEN=16) id
   INTEGER(LONG) saved_seed,inp_seed
   INTEGER clock_rate,clock_1,clock_2
   INTEGER hotcoldread,save,num_upd,p,m,sign,nu,i,hit,num_hit
!!
   INTEGER xyzt
!!

! input variables:
   WRITE (*,'("Lattice size: ",4I5)') NX,NY,NZ,NT
   WRITE (*,'("Enter beta:   ")',ADVANCE='NO')
   READ *,beta
   WRITE (*,'("Enter number of updates:   ")',ADVANCE='NO')
   READ *,num_upd
   WRITE (*,'("Select the starting configuration.  Enter 0 for a &
               &hot start ")')
   WRITE (*,'("1 for a cold start, 2 to read from Disk: ")',ADVANCE='NO')
   READ *,hotcoldread


! other useful variables:
   num_hit=6                  ! number of Metropolis multiple hits
   range_unit=1._REAL8        ! unitary range for the random numbers
   range_small=0.1_REAL8      ! range for the random numbers
   inp_seed=1                 ! input seed for random numbers generator
   zero_matrix%mc=ZERO_m      ! 0 matrix for initializing the staple
   in_filename= 'configuration.in'    
   out_filename='configuration.out'


! initializing system clock
   CALL SYSTEM_CLOCK(clock_1,clock_rate)
   clock_dcl=1._REAL8/clock_rate
   clock_upd=0._REAL8
   clock_plaq=0._REAL8


! initializing random generator and gauge configuration
   SELECT CASE(hotcoldread)
   CASE(0)
     l_seed=.Seed.inp_seed
     DO p=0,1
     DO m=1,4
       ge=.Ggauss.range_unit
       u%uc(p,m)=.Exp.ge
       u%uc(p,m)%parity=p
       u%uc(p,m)%dir=m
     END DO
     END DO
   CASE(1)
     l_seed=.Seed.inp_seed
     unit_matrix%mc=UNIT
     DO p=0,1
     DO m=1,4
       u%uc(p,m)=unit_matrix
       u%uc(p,m)%parity=p
       u%uc(p,m)%dir=m
     END DO
     END DO
   CASE(2)
     CALL read_conf(saved_beta,id,hp,saved_seed,in_filename)    
     IF(inp_seed==0) THEN
       WRITE (*,'("saved_seed=",I15)') saved_seed
       l_seed=.Seed.saved_seed
     ELSE
       l_seed=.Seed.inp_seed
       WRITE (*,'("seed re-initialized")')
     ENDIF
   CASE DEFAULT
     WRITE (*,'("hotcoldread must only be 0,1 or 2")')
     STOP
   END SELECT


   DO i=1,num_upd                        ! Main Loop


! Metropolis update 
      CALL SYSTEM_CLOCK(clock_1)
      DO p=0,1
      DO m=1,4

        ! Staple
        staple=zero_matrix
        staple%parity=p
        staple%dir=m
        DO nu=1,4
          IF(nu.EQ.m) CYCLE
          DO sign=-1,1,2
            staple=staple+((nu*sign).Ushift.u%uc(1-p,m))
          END DO
        END DO

        g_old=u%uc(p,m)
        DO hit=1,num_hit
          plaq_old=g_old.Dot.staple
          ge=.Ggauss.range_small
          ge%parity=p
          ge%dir=m
          g_new=(.Exp.ge)*g_old
          plaq_new=g_new.Dot.staple
          bf_ratio=.Exp.(beta/2._REAL8*(plaq_new-plaq_old))
!          bf_ratio=.Exp.(beta/3._REAL8*(plaq_new-plaq_old))
          rand=.Rand.range_unit
          l_test=rand<bf_ratio
          assign_type='t'; g_old=g_new
!!
!          DO xyzt=0,NXYZT2-1
!             if(context(xyzt)) then
!                g_old%fc(1,1,xyzt)=g_new%fc(1,1,xyzt)
!                g_old%fc(1,2,xyzt)=g_new%fc(1,2,xyzt)
!                g_old%fc(2,1,xyzt)=g_new%fc(2,1,xyzt)
!                g_old%fc(2,2,xyzt)=g_new%fc(2,2,xyzt)
!             end if
!          END DO
!!
        END DO
        u%uc(p,m)=g_old
      END DO
      END DO
      CALL SYSTEM_CLOCK(clock_2)
      clock_upd=clock_upd+(clock_2-clock_1)*clock_dcl


! Plaquette
      CALL SYSTEM_CLOCK(clock_1)
      av_plaq=0._REAL8
      DO p=0,1
      DO m=1,3
        DO nu=m+1,4
          aux=u%uc(p,m).Dot.(nu.Ushift.u%uc(1-p,m))
          av_plaq=av_plaq+aux
        END DO
      END DO
      END DO
      av_plaq=av_plaq/REAL(12*NXYZT,REAL8)
!      av_plaq=av_plaq/REAL(18*NXYZT,REAL8)
      CALL SYSTEM_CLOCK(clock_2)
      clock_plaq=clock_plaq+(clock_2-clock_1)*clock_dcl
     WRITE (*,'("iteration ",I5," av. plaq.= ",F10.6)') i,1._REAL8-av_plaq


   END DO                                ! End Main Loop


! Save configuration on disk
   WRITE (*,'("Save configuration on disk ? (Yes=1, No=0): ")',ADVANCE='NO')
   READ *,save
   IF(save==1) THEN
     WRITE (*,'("saving the configuration")')
     id='conf 0.0.0'
     hp=0.0
     l_seed=.TRUE.
     saved_seed=.Seed.l_seed
     WRITE (*,'("  saved_seed = ",I15)') saved_seed
     CALL write_conf(beta,id,hp,saved_seed,out_filename)
   ENDIF


! Print timing
   WRITE (*,'("Av. upgrade time in microsecs per link",F9.3)') &
            (1000000*clock_upd)/(4*NXYZT*num_upd)
   WRITE (*,'("Av. measure time in microsecs per plaquette",F9.3)') &
            (1000000*clock_plaq)/(6*NXYZT*num_upd)

   END


