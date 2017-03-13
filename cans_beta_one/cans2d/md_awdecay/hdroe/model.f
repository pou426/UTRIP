c======================================================================|
      subroutine model(ro,vx,vy,vz,bx,by,bz,cs2,thini
     &   ,margin,x,ix,y,jx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),vz(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      cs2=1.0d0

      thini=60.d0/180.d0*pi
      rkx=1.d0
      rky=1.d0
      xmax=rkx/cos(thini)
      ymax=rky/sin(thini)
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xmax/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=ix/2+1
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

      jzero=jx/2+1
      y(jzero)=dym(jzero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      ca0=sqrt(10.d0)
      ampaw=sqrt(0.9d0)
      bb0=ca0*sqrt(4.d0*pi)

      do j=1,jx
      do i=1,ix
         s=x(i)*cos(thini)+y(j)*sin(thini)
         phase=2.d0*pi*s
         ro(i,j)  = 1.d0
         vs = 0.0d0
         vn = ampaw*ca0*sin(phase)
         vt = ampaw*ca0*cos(phase)
         vx(i,j)  = vs*cos(thini)-vn*sin(thini)
         vy(i,j)  = vs*sin(thini)+vn*cos(thini)
         vz(i,j)  = vt
         bs = bb0
         bn  = -ampaw*bb0*sin(phase)
         bt  = -ampaw*bb0*cos(phase)
         bx(i,j)  = bs*cos(thini)-bn*sin(thini)
         by(i,j)  = bs*sin(thini)+bn*cos(thini)
         bz(i,j)  = bt
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ca0',ca0)
      call dacputparamd(mf_params,'ampaw',ampaw)
      call dacputparamd(mf_params,'thini',thini)



      return
      end
