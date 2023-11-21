      program test_su2
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c What to Do:
c Monte Carlo Simulation of SU(2) Yang-Mills theory
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calling Subroutines:
c lattice_setup,init,monte
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c I/O: 
c unit=16,file='data.test_su2'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c How to Use:
c let the object read the following data from standard input
c      read(5,*) b, nhit
c      read(5,*) ints,irand0,niter

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc /* global variables and parameters */

      include 'parameters.f'

      common /vmtrx/ u0(nsite,ndim),u1(nsite,ndim),
     c               u2(nsite,ndim),u3(nsite,ndim)
      real*8 u0,u1,u2,u3

      common /vbeta/ b,hitp                                   
      real*8 b,hitp

      common /vdata/ eu
      real*8 eu

      common /viter/ iter                                         
      integer iter

      common /vrand/ iri(250),ir(nvecr)
      integer iri,ir

cc /* local variables */

      integer nhit
      integer ints,irand0
      integer niter

      common /drct/ mup(nsite,ndim),mdn(nsite,ndim),
     c              msite(nx,ny,nz,nt)
      integer mup,mdn,msite

      integer n

cc /* main start */

      read(5,*) b, nhit
      read(5,*) ints,irand0,niter

      open(unit=16,file='data.test_su2',status='unknown')

 1000 format(
     c 1x,'SU(2) Yang-Mills Model'/
     c 1x,'Monte Carlo Simulation'/
     c 1x,'                       '/
     c 1x,'Gauge Configuration   '/
     c 1x,' Lattice               ',i2,' x ',i2,' x ',i2,' x ',i2/
     c 1x,' Beta, Nhit            ',f12.8,1x,i3/
     c 1x,' Initial Condition     ',i1,' (1=cold 2=continued)'/
     c 1x,' Initial Random Number ',i20/
     c 1x,' Number of Iteration   ',i5)

      write(16,1000)   
     c                nx,ny,nz,nt,
     c                b,nhit,
     c                ints,
     c                irand0,
     c                niter
      write(16,*)

      call lattice_setup
      call inits(ints,irand0)

      do n=1,niter
         iter=iter+1 
         call monte(nhit)
         write(6,*) iter, 1.0d0-eu
         write(16,*) iter, 1.0d0-eu
      end do
      
      write(6,*) 
      write(6,*) b, 1.0d0-eu

      stop
      end




