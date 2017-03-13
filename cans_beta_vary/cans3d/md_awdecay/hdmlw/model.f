c======================================================================|
      subroutine model(ro,vx,vy,vz,bx,by,bz,cs2,thini,phini
     &       ,margin,x,ix,y,jx,z,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      cs2=1.d0

      thini=60.d0/180.d0*pi
      phini=30.d0/180.d0*pi
      rkx=1.d0
      rky=1.d0
      rkz=1.d0
      xmax=rkx/(sin(thini)*cos(phini))
      ymax=rky/(sin(thini)*sin(phini))
      zmax=rkz/cos(thini)

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xmax/real(ix-margin*2)
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

      dy0=ymax/real(jx-margin*2)
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

      dz0=zmax/real(kx-margin*2)
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
c     store initial condition into common area
c----------------------------------------------------------------------|
      ca0=sqrt(10.d0)
      ampaw=sqrt(0.9d0)
      bb0=ca0*sqrt(4.d0*pi)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         s=x(i)*sin(thini)*cos(phini)+y(j)*sin(thini)*sin(phini)
     &    +z(k)*cos(thini)
         phase=2.d0*pi*s
         ro(i,j,k)  = 1.d0
         vs = 0.0d0
         vn = ampaw*ca0*sin(phase)
         vt = ampaw*ca0*cos(phase)
         vx(i,j,k)  = vs*sin(thini)*cos(phini)
     &             +vn*(-sin(phini))
     &             +vt*(-cos(thini)*cos(phini))
         vy(i,j,k)  = vs*sin(thini)*sin(phini)
     &             +vn*cos(phini)
     &             +vt*(-cos(thini)*sin(phini))
         vz(i,j,k)  = vs*cos(thini)
     &             +vt*sin(thini)
         bs = bb0
         bn  = -ampaw*bb0*sin(phase)
         bt  = -ampaw*bb0*cos(phase)
         bx(i,j,k)  = bs*sin(thini)*cos(phini)
     &             +bn*(-sin(phini))
     &             +bt*(-cos(thini)*cos(phini))
         by(i,j,k)  = bs*sin(thini)*sin(phini)
     &             +bn*cos(phini)
     &             +bt*(-cos(thini)*sin(phini))
         bz(i,j,k)  = bs*cos(thini)
     &             +bt*sin(thini)

      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ca0',ca0)
      call dacputparamd(mf_params,'ampaw',ampaw)
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'phini',phini)


      return
      end
