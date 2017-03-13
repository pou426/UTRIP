c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm,mu,et,etm
     &       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm
     &       ,margin,ix,jx,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dx(ix),xm(ix),dxm(ix)
      dimension y(jx),dy(jx),ym(jx),dym(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension et(ix,jx,kx),etm(ix,jx,kx)
      double precision mu

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      pi = acos(-1.0d0)
      mu =4.d0*pi

      xmax=5.d0
      ymax=20.d0
      zmax=2.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xmax/dble(max(ix-margin*2,1))
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=min(margin+1,ix)
      x(izero)=dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      call grdrdy2(dx,xm,dxm,x,ix)

      dy0=ymax/dble(max(jx-margin*2,1))
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=min(margin+1,jx)
      y(jzero)=dym(jzero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

      call grdrdy2(dy,ym,dym,y,jx)

      dz0=zmax/dble(max(kx-margin*2,1))
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

c----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c----------------------------------------------------------------------|

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
