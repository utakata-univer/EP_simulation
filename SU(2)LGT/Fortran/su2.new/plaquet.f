      subroutine plaquet(tu0,tu1,tu2,tu3,plaq)

      implicit none
      include 'parameters.f'

      common /drct/ mup(nsite,ndim),mdn(nsite,ndim),
     c              msite(nx,ny,nz,nt)
      integer mup,mdn,msite

      common /tld/ ldir2(ndim1,ndim)
      integer ldir2

      real*8 tu0(nsite,ndim),tu1(nsite,ndim),
     c       tu2(nsite,ndim),tu3(nsite,ndim)
      real*8 plaq

      real*8 wu0(nsite),wu1(nsite),
     c       wu2(nsite),wu3(nsite)

      real*8 wd0(nsite),wd1(nsite),
     c       wd2(nsite),wd3(nsite)

      integer m,mk,ml
      integer k,l,ll

      real*8 dplaq
      integer nplaq

      nplaq=ndim*ndim1*nsite
      dplaq=1.d0/dble(nplaq)

      plaq=0.d0

      do k=1,ndim
         do ll=1,ndim1
            l=ldir2(ll,k)

            do m=1,nsite

               mk=mup(m,k)
               ml=mup(m,l)

               wu0(m) = tu0(m,k)*tu0(mk,l)-tu1(m,k)*tu1(mk,l)            
     c              -tu2(m,k)*tu2(mk,l)-tu3(m,k)*tu3(mk,l)            
               wu1(m) = tu0(m,k)*tu1(mk,l)+tu1(m,k)*tu0(mk,l)            
     c              -tu2(m,k)*tu3(mk,l)+tu3(m,k)*tu2(mk,l)            
               wu2(m) = tu0(m,k)*tu2(mk,l)+tu2(m,k)*tu0(mk,l)            
     c              -tu3(m,k)*tu1(mk,l)+tu1(m,k)*tu3(mk,l)            
               wu3(m) = tu0(m,k)*tu3(mk,l)+tu3(m,k)*tu0(mk,l)            
     c              -tu1(m,k)*tu2(mk,l)+tu2(m,k)*tu1(mk,l)            


               wd0(m) = tu0(m,l)*tu0(ml,k)-tu1(m,l)*tu1(ml,k)
     c              -tu2(m,l)*tu2(ml,k)-tu3(m,l)*tu3(ml,k)            
               wd1(m) =tu0(m,l)*tu1(ml,k)+tu1(m,l)*tu0(ml,k)            
     c              -tu2(m,l)*tu3(ml,k)+tu3(m,l)*tu2(ml,k)            
               wd2(m) =tu0(m,l)*tu2(ml,k)+tu2(m,l)*tu0(ml,k)            
     c              -tu3(m,l)*tu1(ml,k)+tu1(m,l)*tu3(ml,k)            
               wd3(m) =tu0(m,l)*tu3(ml,k)+tu3(m,l)*tu0(ml,k)            
     c              -tu1(m,l)*tu2(ml,k)+tu2(m,l)*tu1(ml,k)            
                                                                        
            end do

            do m=1,nsite
             plaq=plaq
     c           +wu0(m)*wd0(m)+wu1(m)*wd1(m)
     c           +wu2(m)*wd2(m)+wu3(m)*wd3(m)
            end do

         end do
      end do

      plaq=1.d0-plaq*dplaq

      return
      end


