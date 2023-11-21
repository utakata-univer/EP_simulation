      subroutine inits(ints,irand0)

      implicit none
      include 'parameters.f'

      common /viter/ iter
      integer iter

      common /vrand/ iri(250),ir(nvecr) 
      integer iri,ir

      common /vmtrx/ u0(nsite,ndim),u1(nsite,ndim),
     c               u2(nsite,ndim),u3(nsite,ndim)
      real*8 u0,u1,u2,u3

      integer ints,irand0
      integer m,l,nran

      if (ints.eq.1) then
         iter=0

         iri(1) = irand0
         do m = 2,250
            iri(m) = iand(nzz,nr*iri(m-1))
         end do

         do nran = 1,1000
            do m = 1,250
               iri(m) = iand(nzz,nrm*iri(m))
            end do
         end do

         do l = 1,ndim
            do m = 1,nsite

               u0(m,l) = 1.0
               u1(m,l) = 0.0
               u2(m,l) = 0.0
               u3(m,l) = 0.0

            end do
         end do

      else
         if (ints.eq.2) then

            open(unit=10,file='config.data',status='unknown')
            open(unit=12,file='rnd.data',status='unknown')

            read(10,*) u0,u1,u2,u3
            read(12,*) iter,iri,ir

            close(10)
            close(12)

         endif
      endif

      return
      end
                                                                       
