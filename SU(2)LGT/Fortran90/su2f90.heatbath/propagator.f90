!  Program Qcdf90_propagator version 4.0.1

!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.


PROGRAM Qcdf90_propagator

   USE precisions
   USE constants
   USE global_module
   USE field_algebra
   USE generator_algebra
   USE conditionals
   USE shift
   USE dirac
   USE random_numbers
   USE assign_mixed
   USE assign_isotype1
   USE assign_isotype2


   IMPLICIT NONE
   TYPE(fermi_field):: psi,chi,grad,h,m_h,mp_m_h
   REAL(REAL8) clock_dcl,clock_cg
   REAL(REAL8) kappa,tolerance,residue,saved_beta,hp
   REAL(REAL8) alpha,old_alpha,g_2,g_old_2,beta_cg,old_beta_cg
   REAL(REAL8) h_a_h,norm_psi
   CHARACTER(LEN=64) in_filename
   CHARACTER(LEN=16) id
   INTEGER(LONG) saved_seed
   INTEGER clock_rate,clock_1,clock_2
   INTEGER iter,nsteps,niter,init_niter,stop_flag,init_stop_flag
   INTEGER i,xyzt,s

! input variables:
   WRITE (*,'("Enter kappa:   ")',ADVANCE='NO')
   READ *,kappa
   WRITE (*,'("Enter max numbers of cg steps:   ")',ADVANCE='NO')
   READ *,nsteps
   WRITE (*,'("Enter tolerance:   ")',ADVANCE='NO')
   READ *,tolerance
                              ! the conjugated gradient will
                              ! run until the residue<tolerance
                              ! or for a maximum of nsteps


! other useful variables:
   in_filename= 'configuration.in'
   init_stop_flag=2
   init_niter=4


! gauge configuration is read from the disk
   CALL read_conf(saved_beta,id,hp,saved_seed,in_filename)


! the source chi (in the even sites) is set arbitrarly in this example
   DO i=1,3
   DO s=1,4
   DO xyzt=0,NXYZT2-1
     chi%fc(i,xyzt,s)=(1._REAL8,1._REAL8)/SQRT(REAL(24*NXYZT2,REAL8))
   END DO
   END DO
   END DO
   chi%parity=0


! psi must be initialized as the starting trial configuration.
! the simpliest choice is psi=chi
   psi=chi


! initializing system clock
   CALL SYSTEM_CLOCK(clock_1,clock_rate)
   clock_dcl=1._REAL8/clock_rate


!  Calculate psi as the solution of: M*psi=chi,
!  where M is the fermion matrix, and chi is a given source.
!  The residue is printed to monitor the convergence.

    stop_flag=init_stop_flag
    niter=init_niter
    iter=0
    m_h=psi-(kappa**2)*(.Dirac.(.Dirac.psi))
    mp_m_h=m_h-(kappa**2)*(.Xdirac.(.Xdirac.m_h))
    grad=chi-mp_m_h
    g_2=grad*grad
    h=grad
    norm_psi=psi*psi
    residue=g_2/norm_psi
    WRITE (*,'("residue= ",F20.16," at step:",I5)') residue,iter
    old_alpha=0._REAL8
    old_beta_cg=1._REAL8

    DO iter=1,nsteps
      m_h=h-(kappa**2)*(.Dirac.(.Dirac.h))
      h_a_h=m_h*m_h
      beta_cg=g_2/h_a_h
      psi=psi+beta_cg*h
      norm_psi=psi*psi
      IF(mod(iter,niter)==0 .AND. g_2/norm_psi<tolerance) THEN
        stop_flag=stop_flag-1
        m_h=psi-(kappa**2)*(.Dirac.(.Dirac.psi))
        mp_m_h=m_h-(kappa**2)*(.Xdirac.(.Xdirac.m_h))
        grad=chi-mp_m_h
        g_2=grad*grad
        h=grad
        g_old_2=g_2
        residue=g_2/norm_psi
        WRITE (*,'("residue= ",F20.16," at step:",I5)') residue,iter
        IF(stop_flag == 0) EXIT
      ELSE
        mp_m_h=m_h-(kappa**2)*(.Xdirac.(.Xdirac.m_h))
        grad=grad-beta_cg*mp_m_h
        g_old_2=g_2
        g_2=grad*grad
        alpha=g_2/g_old_2
        h=grad+alpha*h
        norm_psi=psi*psi
        residue=g_2/norm_psi
        IF(mod(iter,niter) == 0) THEN
           WRITE (*,'("residue= ",F20.16," at step:",I5)') residue,iter
        ENDIF
        old_beta_cg=beta_cg
        old_alpha=alpha
      END IF

    END DO
    CALL SYSTEM_CLOCK(clock_2)
    clock_cg=(clock_2-clock_1)*clock_dcl


!test solution:
   m_h=psi-(kappa**2)*(.Dirac.(.Dirac.psi))
   mp_m_h=m_h-(kappa**2)*(.Xdirac.(.Xdirac.m_h))
   grad=chi-mp_m_h
   norm_psi=psi*psi
   g_2=grad*grad
   residue=g_2/norm_psi
   WRITE (*,'("final residue= ",F20.16)') residue


! Print timing
   WRITE (*,'("Cg time per iteration per link in microsecs",F9.3)') &
         (1000000*clock_cg)/(iter*4*NXYZT)


   END


