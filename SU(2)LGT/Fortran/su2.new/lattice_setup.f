      subroutine lattice_setup
      implicit none

      include 'parameters.f'

      common /drct/ mup(nsite,ndim),mdn(nsite,ndim),
     c              msite(nx,ny,nz,nt)
      integer mup,mdn,msite

      common /site/ ix(nsite),iy(nsite),iz(nsite),it(nsite)
      integer ix,iy,iz,it

      common /oe/ moeup(nvect,ndim,2),moedn(nvect,ndim,2),
     c            ms(nvect,2)
      integer moeup,moedn,ms

      common /tld/ ldir2(ndim1,ndim)
      integer ldir2

      common /vtwst/ zt(nsite,ndim,ndim)                                
      real*8 zt

      integer mu(nt,ndim),md(nt,ndim),mcycle(0:n0+1),isize(ndim)
      integer ipower(ndim)

      integer i,j
      integer m,n,i1,i2,i3,i4
      integer m4,mm4,mupl,mdnl
      integer m3,mm3,mupk,mdnk
      integer m2,mm2,mupj,mdnj
      integer m1,mm1,mupi,mdni
      integer ioe

c
      ldir2(1,1)  =  2
      ldir2(2,1)  =  3
      ldir2(3,1)  =  4
      ldir2(1,2)  =  3
      ldir2(2,2)  =  4
      ldir2(3,2)  =  1
      ldir2(1,3)  =  4
      ldir2(2,3)  =  1
      ldir2(3,3)  =  2
      ldir2(1,4)  =  1
      ldir2(2,4)  =  2
      ldir2(3,4)  =  3

      isize(1)  =  nx
      isize(2)  =  ny
      isize(3)  =  nz
      isize(4)  =  nt

      ipower(1) =  1
      ipower(2) =  nx
      ipower(3) =  nx*ny
      ipower(4) =  nx*ny*nz

      do 1 n  = 1,ndim
      do 1 m  = 1,nt
      mu(m,n)  =  mod(m,isize(n))+1
  1   md(m,n)  =  mod(m-2+isize(n),isize(n))+1

      do 11 m = 0,n0+1
  11  mcycle(m) = mod(m-1+n0,n0) + 1

      do 2 i4 = 1,nt

      m4  = (i4-1)*ipower(4)
      mm4 = (i4-1)*ipower(4)/2
      mupl = (mu(i4,4)-i4)*ipower(4)/2
      mdnl = (md(i4,4)-i4)*ipower(4)/2

      do 2 i3 = 1,nz

      m3  = m4  + (i3-1)*ipower(3)
      mm3 = mm4 + (i3-1)*ipower(3)/2
      mupk = (mu(i3,3)-i3)*ipower(3)/2
      mdnk = (md(i3,3)-i3)*ipower(3)/2

      do 2 i2 = 1,ny

      m2  = m3  + (i2-1)*ipower(2)
      mm2 = mm3 + (i2-1)*ipower(2)/2
      mupj = (mu(i2,2)-i2)*ipower(2)/2
      mdnj = (md(i2,2)-i2)*ipower(2)/2

      ioe = mod(i2+i3+i4,2)

      do 2 i1 = 1,n0

      m1  = m2 +2*i1
      mm1 = mm2+i1

      moeup(mm1,1,1) = mm1 + mcycle(i1+1-ioe)-i1
      moedn(mm1,1,1) = mm1 + mcycle(i1  -ioe)-i1
      moeup(mm1,1,2) = mm1 + mcycle(i1+ioe  )-i1
      moedn(mm1,1,2) = mm1 + mcycle(i1+ioe-1)-i1

      moeup(mm1,2,1) = mm1 + mupj
      moedn(mm1,2,1) = mm1 + mdnj
      moeup(mm1,2,2) = mm1 + mupj
      moedn(mm1,2,2) = mm1 + mdnj

      moeup(mm1,3,1) = mm1 + mupk
      moedn(mm1,3,1) = mm1 + mdnk
      moeup(mm1,3,2) = mm1 + mupk
      moedn(mm1,3,2) = mm1 + mdnk

      moeup(mm1,4,1) = mm1 + mupl
      moedn(mm1,4,1) = mm1 + mdnl
      moeup(mm1,4,2) = mm1 + mupl
      moedn(mm1,4,2) = mm1 + mdnl

      ms(mm1,1) = m1-ioe
      ms(mm1,2) = m1+ioe-1

  2   continue

      do 4 i4 = 1,nt

      m4   = (i4-1)*ipower(4)
      mupl = (mu(i4,4)-i4)*ipower(4)
      mdnl = (md(i4,4)-i4)*ipower(4)

      do 4 i3 = 1,nz

      m3   = m4 + (i3-1)*ipower(3)
      mupk = (mu(i3,3)-i3)*ipower(3)
      mdnk = (md(i3,3)-i3)*ipower(3)

      do 4 i2 = 1,ny

      m2   = m3 + (i2-1)*ipower(2)
      mupj = (mu(i2,2)-i2)*ipower(2)
      mdnj = (md(i2,2)-i2)*ipower(2)

      do 4 i1 = 1,nx

      m1   = m2 + i1
      mupi = mu(i1,1)-i1
      mdni = md(i1,1)-i1

      mup(m1,1) = m1+mupi
      mup(m1,2) = m1+mupj
      mup(m1,3) = m1+mupk
      mup(m1,4) = m1+mupl
      mdn(m1,1) = m1+mdni
      mdn(m1,2) = m1+mdnj
      mdn(m1,3) = m1+mdnk
      mdn(m1,4) = m1+mdnl

      msite(i1,i2,i3,i4)=m1
      ix(m1)=i1
      iy(m1)=i2
      iz(m1)=i3
      it(m1)=i4

  4   continue

      do 5 i4 = 1,nt                                                    

         m4 = (i4-1)*ipower(4)                                             

         do 5 i3 = 1,nz                                                    

            m3 = m4 + (i3-1)*ipower(3)                                        

            do 5 i2 = 1,ny                                                    

               m2 = m3 + (i2-1)*ipower(2)

               do 5 i1 = 1,nx

                  m1 = m2 + i1

                  do 6 i = 1,ndim 
                     zt(m1,i,i) = 0.d0  
 6                continue
                  do 7 j = 1,ndim
                     do 7 i = 1,ndim
                        if (i.eq.j) goto 7 
                        zt(m1,i,j) = 1.d0 
 7                   continue

  5   continue                                           

      return
      end

