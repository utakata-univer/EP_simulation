      subroutine rnd

      implicit none
      include 'parameters.f'

      common /vrand/ iri(250),ir(nvecr) 
      integer iri,ir
      common /srand/ iy(128)
      integer iy

      integer i,j

      do 20 j = 0,nvect/128-1

      do 21 i = 1,128
      iy(i) = ieor(iri(i),iri(i+103))
      ir(j*128+i) = iy(i)
 21   continue

      do 22 i = 1,122
 22   iri(i) = iri(i+128)

      do 23 i = 1,128
 23   iri(i+122) = iy(i)

 20   continue

      return
      end
