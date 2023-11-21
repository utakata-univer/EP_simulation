!  Program Qcdf90_quenched version 4.0.0

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  May 1996
!  This program may be freely copied and used as long as this notice
!  is retained.

! su2 version  2006.08.28 by y.k.
! heatbath version
PROGRAM su2f90_heatbath

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

   TYPE(real_field):: a0,a1,a2,a3,a_norm,cos_theta,sin_theta,phi
   TYPE(real_field):: norm_staple,bk,norm_g
   TYPE(real_field):: a0_tmp,z,a0_old
   TYPE(gauge_field):: g_unit,g_sigma1,g_sigma2,g_sigma3
   TYPE(matrix):: sigma1_matrix,sigma2_matrix,sigma3_matrix
   INTEGER xyzt,q,n

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

   unit_matrix%mc=UNIT
   sigma1_matrix%mc=SIGMA1
   sigma2_matrix%mc=SIGMA2
   sigma3_matrix%mc=SIGMA3

   g_unit = unit_matrix
   g_sigma1 = sigma1_matrix
   g_sigma2 = sigma2_matrix
   g_sigma3 = sigma3_matrix
   

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

! Heat Bath update 
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
        norm_staple=.Sqrt.((staple.Dot.staple)/2._REAL8)
!
        bk=beta*norm_staple

        g_old =(.Adj.(staple/norm_staple))*u%uc(p,m)
        a0_old=(.Tr.g_old)/2._REAL8

        a0=a0_old
        DO hit=1,num_hit
!         z=(.Exp.(bk)-.Exp.(.Minus.bk))*.Rand.range_unit+.Exp.(.Minus.bk)
!         a0_tmp=(.Log.(z))/bk
         z=(.Exp.(2._REAL8*bk)-1._REAL8)*.Rand.range_unit+1._REAL8
         a0_tmp=(.Log.(z))/bk - 1._REAL8
         a0_tmp%parity=p

         a_norm=.Sqrt.(1._REAL8-a0_tmp*a0_tmp)
         rand=.Rand.range_unit
         l_test=rand<a_norm
!         assign_type='t'; a0=a0_tmp
         DO xyzt=0,NXYZT2-1
            if(context(xyzt)) then
               a0%fc(xyzt)=a0_tmp%fc(xyzt)
            end if
          END DO
        END DO

        a_norm=.Sqrt.(1._REAL8-a0*a0)
        cos_theta=2._REAL8*(.Rand.range_unit)-1._REAL8
        sin_theta=.Sqrt.(1._REAL8-cos_theta*cos_theta)
        phi=TWOPI*.Rand.range_unit

        a3=a_norm*cos_theta 
        a1=a_norm*sin_theta*.Cos.phi
        a2=a_norm*sin_theta*.Sin.phi

        g_new=a0*g_unit+.I.(a1*g_sigma1+a2*g_sigma2+a3*g_sigma3)
        g_new%parity=p
        g_new%dir=m

        u%uc(p,m)=(staple/norm_staple)*g_new

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


