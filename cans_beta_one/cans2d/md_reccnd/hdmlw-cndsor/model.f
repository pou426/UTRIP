c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm,rkap0
     &           ,et,etm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)
      dimension et(ix,jx),etm(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      rkap0=1.0d0
      gm=5.d0/3.d0
      pi = acos(-1.0d0)
      xmax=5.d0
      ymax=20.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xmax/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      x(izero)=dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=ymax/dble(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin+1
      y(jzero)=dym(jzero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c   resistivity
c----------------------------------------------------------------------|
      etai=0.1d0
c     etai=0.0
      retai=0.8d0

      do j=1,jx
      do i=1,ix
         rr=sqrt(x(i)**2+y(j)**2)
         if (rr.lt.retai) then
           et(i,j)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
         rr=sqrt((x(i)+dxm(i)/2)**2+(y(j)+dym(j)/2)**2)
         if (rr.lt.retai) then
           etm(i,j)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
      enddo
      enddo
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      betai=10.d0

      do j=1,jx
      do i=1,ix
         bx(i,j)=0.d0
         by(i,j)=sqrt(betai*pi*8/gm)*tanh(x(i)/0.5)
         bz(i,j)=sqrt(betai*pi*8/gm)/cosh(x(i)/0.5)
         pr(i,j)=1/gm
         ro(i,j)=pr(i,j)*gm
         vx(i,j)=0.d0
         vy(i,j)=0.d0
         vz(i,j)=0.d0
      enddo
      enddo


c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'betai',betai)
      call dacputparamd(mf_params,'rkap0',rkap0)


      return
      end
