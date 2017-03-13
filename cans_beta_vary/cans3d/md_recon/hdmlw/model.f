c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &           ,et,etm,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension et(ix,jx,kx),etm(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      pi = acos(-1.0d0)

      xmax=5.d0
      ymax=20.d0
      zmax=2.d0
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

      dz0=zmax/dble(kx-margin*2)
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=margin+1
      z(kzero)=dzm(kzero)/2.d0
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo

c----------------------------------------------------------------------|
c   resistivity
c----------------------------------------------------------------------|
      etai=0.03d0
      retai=0.8d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
c        rr=sqrt(x(i)**2+y(j)**2+z(k)**2)
         rr=sqrt(x(i)**2+y(j)**2)
         if (rr.lt.retai) then
           et(i,j,k)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
c      rr=sqrt((x(i)+dxm(i)/2)**2+(y(j)+dym(j)/2)**2+(z(k)+dzm(k)/2)**2)
       rr=sqrt((x(i)+dxm(i)/2)**2+(y(j)+dym(j)/2)**2)
         if (rr.lt.retai) then
           etm(i,j,k)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
      enddo
      enddo
      enddo
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      betai=10.d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         bx(i,j,k)=0.d0
         by(i,j,k)=sqrt(betai*pi*8/gm)*tanh(x(i)/0.5d0)
         bz(i,j,k)=0.d0
         pr(i,j,k)=1/gm*(1+betai)-by(i,j,k)**2/pi/8.d0
         ro(i,j,k)=pr(i,j,k)*gm

         vx(i,j,k)=0.0d0
         vy(i,j,k)=0.0d0
         vz(i,j,k)=0.0d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'betai',betai)


      return
      end
