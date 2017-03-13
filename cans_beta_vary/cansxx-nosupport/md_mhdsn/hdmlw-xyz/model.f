c======================================================================|
      subroutine model(ro,pr,vr,vph,vz,br,bph,bz,gm,mu
     &   ,r,dr,rm,drm,ph,dph,phm,dphm,z,dz,zm,dzm
     &   ,margin,ix,jx,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension r(ix),dr(ix),rm(ix),drm(ix)
      dimension ph(jx),dph(jx),phm(jx),dphm(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vr(ix,jx,kx),vph(ix,jx,kx),vz(ix,jx,kx)
      dimension br(ix,jx,kx),bph(ix,jx,kx),bz(ix,jx,kx)
      double precision mu

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      pi=acos(-1.0d0)
      mu =4.d0*pi

      rmin=0.04d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/dble(max(ix-margin*2,1))
      do i=1,ix
         drm(i)=dx0
      enddo
       
      izero=min(margin+1,ix)
      r(izero)=rmin+drm(izero)/2.d0
      do i=izero+1,ix
         r(i) = r(i-1)+drm(i-1)
      enddo
      do i=izero-1,1,-1
         r(i) = r(i+1)-drm(i)
      enddo

      call grdrdy2(dr,rm,drm,r,ix)

      dy0=2.d0*pi/dble(max(jx-margin*2,1))
      do j=1,jx
         dphm(j)=dy0
      enddo
       
      jzero=min(margin+1,jx)
      ph(jzero)=dphm(jzero)/2.d0
      do j=jzero+1,jx
         ph(j) = ph(j-1)+dphm(j-1)
      enddo
      do j=jzero-1,1,-1
         ph(j) = ph(j+1)-dphm(j)
      enddo

      call grdrdy2(dph,phm,dphm,ph,jx)

      dz0=1.d0/dble(max(kx-margin*2,1))
      do k=1,kx
         dzm(k)=dz0
      enddo
       
      kzero=min(margin+1,kx)
      z(kzero)=dzm(kzero)/2.d0
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo

      call grdrdy2(dz,zm,dzm,z,kx)

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      prism=1.d-8
      wexp=0.1d0

      betai=1.0d6
      b0=sqrt(2.d0*mu*prism*betai)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k) = 1.0d0
         vr(i,j,k) = 0.0d0
         vph(i,j,k) = 0.0d0
         vz(i,j,k) = 0.0d0
         ss=sqrt(r(i)**2+z(k)**2)
         pr(i,j,k) = prism+(1.d0/gm-prism)*exp(-(ss/wexp)**2)
         br(i,j,k) = 0.0d0
         bph(i,j,k) = 0.0d0
         bz(i,j,k) = b0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'wexp',wexp)
      call dacputparamd(mf_params,'prism',prism)
      call dacputparamd(mf_params,'betai',betai)


      
      return
      end
