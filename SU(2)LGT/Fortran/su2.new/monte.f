      subroutine monte(nhit)                                      
                                                                  
      implicit none
      include 'parameters.f'
                                                                        
      common /vrand/ iri(250),ir(nvecr) 
      integer iri,ir
                                                                  
      common /vmtrx/ u0(nsite,ndim),u1(nsite,ndim),
     c               u2(nsite,ndim),u3(nsite,ndim)
      real*8 u0,u1,u2,u3

      common /vbeta/ b,hitp                                   
      real*8 b,hitp

      common /vdata/ eu
      real*8 eu
                                                                  
      common /viter/ iter                                         
      integer iter

      common /drct/ mup(nsite,ndim),mdn(nsite,ndim),
     c              msite(nx,ny,nz,nt)
      integer mup,mdn,msite

      common /oe/ moeup(nvect,ndim,2),moedn(nvect,ndim,2),
     c            ms(nvect,2)
      integer moeup,moedn,ms

      common /tld/ ldir2(ndim1,ndim)
      integer ldir2

      common /vtwst/ zt(nsite,ndim,ndim)                                
      real*8 zt

      common /smonr/ w10(nvect,ndim,ndim),w11(nvect,ndim,ndim),          
     c               w12(nvect,ndim,ndim),w13(nvect,ndim,ndim),          
     c               w20(nvect,ndim,ndim),w21(nvect,ndim,ndim),          
     c               w22(nvect,ndim,ndim),w23(nvect,ndim,ndim),          
     c               w0(nvect),w1(nvect),w2(nvect),w3(nvect),            
     c               wab(nvect), bred(nvect), bexp(nvect), d0(nvect),    
     c               dt(nvect),  d3(nvect),   rad2(nvect)               
      real*8 w10,w11,w12,w13,w20,w21,w22,w23
      real*8 w0,w1,w2,w3
      real*8 wab,bred,bexp,d0,dt,d3,rad2
                                                                        
      common /smoni/ idf(nvect)                                         
      integer idf

      integer nhit

      integer nnv,nnv2,l1,l2,ll
      integer m,mi,m2,m4,m6,m5
      integer nht
      real*8 ztt,ztc
      real*8 rad
      real*8 teta,d1,d2

      hitp = 0.d0
      eu = 0.d0
                                                                        
      do 51 nnv = 1,2

         nnv2 = mod(nnv,2) + 1
                                                                        
         do 101 l2 = 2,ndim                                                
            do 101 l1 = 1,l2-1                                                
                                                                        
               do 102 m = 1,nvect
                  mi = ms(m,nnv)
                  m2 = mup(mi,l1) 
                  m4 = mup(mi,l2)
                  m6 = mdn(mi,l2)
                                                                        
      w10(m,l1,l2) = u0(m2,l2)*u0(m4,l1)+u1(m2,l2)*u1(m4,l1)            
     c              +u2(m2,l2)*u2(m4,l1)+u3(m2,l2)*u3(m4,l1)            
                                                                        
      w11(m,l1,l2) =-u0(m2,l2)*u1(m4,l1)+u1(m2,l2)*u0(m4,l1)            
     c              +u2(m2,l2)*u3(m4,l1)-u3(m2,l2)*u2(m4,l1)            
                                                                        
      w12(m,l1,l2) =-u0(m2,l2)*u2(m4,l1)-u1(m2,l2)*u3(m4,l1)            
     c              +u2(m2,l2)*u0(m4,l1)+u3(m2,l2)*u1(m4,l1)            
                                                                        
      w13(m,l1,l2) =-u0(m2,l2)*u3(m4,l1)+u1(m2,l2)*u2(m4,l1)            
     c              -u2(m2,l2)*u1(m4,l1)+u3(m2,l2)*u0(m4,l1)            
                                                                        
      w20(m,l1,l2) = u0(m6,l1)*u0(m6,l2)+u1(m6,l1)*u1(m6,l2)            
     c              +u2(m6,l1)*u2(m6,l2)+u3(m6,l1)*u3(m6,l2)            
                                                                        
      w21(m,l1,l2) = u0(m6,l1)*u1(m6,l2)-u1(m6,l1)*u0(m6,l2)            
     c              +u2(m6,l1)*u3(m6,l2)-u3(m6,l1)*u2(m6,l2)            
                                                                        
      w22(m,l1,l2) = u0(m6,l1)*u2(m6,l2)-u1(m6,l1)*u3(m6,l2)            
     c              -u2(m6,l1)*u0(m6,l2)+u3(m6,l1)*u1(m6,l2)            
                                                                        
      w23(m,l1,l2) = u0(m6,l1)*u3(m6,l2)+u1(m6,l1)*u2(m6,l2)            
     c              -u2(m6,l1)*u1(m6,l2)-u3(m6,l1)*u0(m6,l2)            
                                                                        
  102 continue                                                          
  101 continue                                                          
                                                                        
      do 103 l2 = 1,ndim1                                               
      do 103 l1 = l2+1,ndim                                             
                                                                        
c*voption vec                                                            
*vocl loop,novrec
      do 104 m = 1,nvect                                                
                                                                        
      m5 = moedn(moeup(m,l1,nnv),l2,nnv2)                               
                                                                        
      w10(m,l1,l2) = w10(m,l2,l1)                                       
      w11(m,l1,l2) =-w11(m,l2,l1)                                       
      w12(m,l1,l2) =-w12(m,l2,l1)                                       
      w13(m,l1,l2) =-w13(m,l2,l1)                                       
      w20(m,l1,l2) = w20(m5,l2,l1)                                      
      w21(m,l1,l2) =-w21(m5,l2,l1)                                      
      w22(m,l1,l2) =-w22(m5,l2,l1)                                      
      w23(m,l1,l2) =-w23(m5,l2,l1)                                      
                                                                        
  104 continue                                                          
  103 continue                                                          
                                                                        
      do 151 l1 = 1,ndim                                                
                                                                        
      do 106 m = 1,nvect                                                
                                                                        
      w0(m) = 0.                                                        
      w1(m) = 0.                                                        
      w2(m) = 0.                                                        
  106 w3(m) = 0.                                                        
                                                                        
      do 107 ll = 1,ndim1                                               
                                                                        
      l2 = ldir2(ll,l1)                                                 
                                                                        
      do 108 m = 1,nvect                                                
                                                                        
      mi = ms(m,nnv)                                                    
                                                                        
      m5 = mdn(mup(mi,l1),l2)                                           
      m6 = mdn(mi,l2)                                                   
                                                                        
      ztt = zt(mi,l1,l2)                                                
      ztc = zt(m6,l2,l1)                                                
                                                                        
      w0(m) = w0(m)                                                     
     c      + ztt*(+w10(m,l1,l2)*u0(mi,l2)+w11(m,l1,l2)*u1(mi,l2)       
     c             +w12(m,l1,l2)*u2(mi,l2)+w13(m,l1,l2)*u3(mi,l2))      
     c      + ztc*(+u0(m5,l2)*w20(m,l1,l2)+u1(m5,l2)*w21(m,l1,l2)       
     c             +u2(m5,l2)*w22(m,l1,l2)+u3(m5,l2)*w23(m,l1,l2))      
                                                                        
      w1(m) = w1(m)                                                     
     c      + ztt*(+w10(m,l1,l2)*u1(mi,l2)-w11(m,l1,l2)*u0(mi,l2)       
     c             -w12(m,l1,l2)*u3(mi,l2)+w13(m,l1,l2)*u2(mi,l2))      
     c      + ztc*(-u0(m5,l2)*w21(m,l1,l2)+u1(m5,l2)*w20(m,l1,l2)       
     c             -u2(m5,l2)*w23(m,l1,l2)+u3(m5,l2)*w22(m,l1,l2))      
                                                                        
      w2(m) = w2(m)                                                     
     c      + ztt*(+w10(m,l1,l2)*u2(mi,l2)+w11(m,l1,l2)*u3(mi,l2)       
     c             -w12(m,l1,l2)*u0(mi,l2)-w13(m,l1,l2)*u1(mi,l2))      
     c      + ztc*(-u0(m5,l2)*w22(m,l1,l2)+u1(m5,l2)*w23(m,l1,l2)       
     c             +u2(m5,l2)*w20(m,l1,l2)-u3(m5,l2)*w21(m,l1,l2))      
                                                                        
      w3(m) = w3(m)                                                     
     c      + ztt*(+w10(m,l1,l2)*u3(mi,l2)-w11(m,l1,l2)*u2(mi,l2)       
     c             +w12(m,l1,l2)*u1(mi,l2)-w13(m,l1,l2)*u0(mi,l2))      
     c      + ztc*(-u0(m5,l2)*w23(m,l1,l2)-u1(m5,l2)*w22(m,l1,l2)       
     c             +u2(m5,l2)*w21(m,l1,l2)+u3(m5,l2)*w20(m,l1,l2))      
                                                                        
  108 continue                                                          
  107 continue                                                          
                                                                        
      do 111 m = 1,nvect                                                
                                                                        
      wab(m) = dsqrt(w0(m)**2+w1(m)**2+w2(m)**2+w3(m)**2)                
                                                                        
      bred(m) = 2.d0*b*wab(m)                                           
      bexp(m) = 1.d0-dexp(-bred(m))                                      
                                                                        
      wab(m) = 1.d0/wab(m)                                              
                                                                        
      idf(m) = 0                                                        
                                                                        
  111 continue                                                          
                                                                        
      do 112 nht = 1,nhit                                               
                                                                        
      call rnd                                                          
      do 113 m = 1,nvect                                                
      if (idf(m).eq.1) goto 113                                         
      dt(m) = 1.d0+2.d0*dlog(1.d0-azz*ir(m)*bexp(m))/bred(m)            
  113 continue                                                          
                                                                        
      call rnd                                                          
      do 114 m = 1,nvect                                                
      if (idf(m).eq.1) goto 114                                         
      if ((azz*ir(m))**2.gt.1.d0-dt(m)**2) goto 114                     
                                                                        
      idf(m) = 1                                                        
      d0(m) = dt(m)                                                     
                                                                        
  114 continue                                                          
  112 continue                                                          
                                                                        
      call rnd                                                          
      do 115 m = 1,nvect                                                
      if (idf(m).eq.0) goto 115                                         
      hitp = hitp + 1.d0                                                
      rad = 1.d0 - d0(m)**2                                             
      d3(m) = dsqrt(dabs(rad))*(2.d0*azz*ir(m)-1.d0)                      
      rad2(m) = dsqrt(dabs(rad-d3(m)**2))                                 
 115  continue                                                          
                                                                        
      call rnd                                                          
c*voption vec                                                            
*vocl loop,novrec
      do 116 m = 1,nvect                                                
                                                                        
      mi = ms(m,nnv)                                                    
                                                                        
      if (idf(m).eq.0) goto 117                                         
                                                                        
      teta = pi2*azz*ir(m)                                              
      d1 = rad2(m)*dcos(teta)                                            
      d2 = rad2(m)*dsin(teta)                                            
                                                                        
      u0(mi,l1) = wab(m)*(d0(m)*w0(m)-d1*w1(m)-d2*w2(m)-d3(m)*w3(m))    
      u1(mi,l1) = wab(m)*(d0(m)*w1(m)+d1*w0(m)-d2*w3(m)+d3(m)*w2(m))    
      u2(mi,l1) = wab(m)*(d0(m)*w2(m)+d1*w3(m)+d2*w0(m)-d3(m)*w1(m))    
      u3(mi,l1) = wab(m)*(d0(m)*w3(m)-d1*w2(m)+d2*w1(m)+d3(m)*w0(m))    
                                                                        
  117 continue                                                          
                                                                        
      eu = eu + u0(mi,l1)*w0(m) + u1(mi,l1)*w1(m)                       
     c        + u2(mi,l1)*w2(m) + u3(mi,l1)*w3(m)                       
                                                                        
  116 continue                                                          
                                                                        
  151 continue                                                          
  51  continue                                                          
                                                                        
      eu = eu/dble(6*nlink)                                            
      hitp = hitp/dble(nlink)                                          
                                                                        
      return                                                            
      end                                                               

