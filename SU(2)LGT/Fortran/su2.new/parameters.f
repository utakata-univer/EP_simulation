      integer ns,nx,ny,nz,nt
      parameter  (ns=8)
      parameter  (nx=ns, ny=ns, nz=ns, nt=8)
c      parameter  (nx=ns, ny=ns, nz=ns, nt=ns*2)
      integer n0
      parameter  (n0=nx/2)
      integer ndim
      parameter  (ndim=4)
      integer ndim1
      parameter  (ndim1=ndim-1)
      integer nsite
      parameter  (nsite=nx*ny*nz*nt)
      integer nlink
      parameter (nlink=nsite*ndim)
      integer nvect
      parameter  (nvect=nsite/2)


      integer nvecr
      parameter  (nvecr=nvect)

c      integer ipp,iqq,ipq
c      parameter  (ipp=607, iqq=273, ipq=ipp-iqq)
c      integer nvecr
c      parameter  (nvecr=(nvect/ipq+1)*ipq)
c      integer nsitr
c      parameter  (nsitr=(nsite/ipq+1)*ipq)
      
      integer nr,nrm,nzz
      parameter(nr = 48828125,nrm = nr**(nsite+1),nzz = 2147483647)
      real*8 azz
      parameter(azz = 1.d0/2147483647.d0)

      real*8 pi2
      parameter(pi2=6.28318530717959)
      real*8 dpi2
      parameter(dpi2=1.0d0/pi2)

      

