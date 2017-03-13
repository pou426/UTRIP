c======================================================================|
      subroutine model(ro,vx,vy,bx,by,cs2,g0,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      cs2=1.0
      g0=1.0

      rlambda=10.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      xmax=rlambda

      dx0=xmax/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      x(izero)=0.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=xmax/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin+1
      y(jzero)=0.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      amp=0.1d0
      betai=0.1d0

      do j=1,jx
      do i=1,ix
         ro(i,j)=1.0+amp*sin(2*pi*x(i)/rlambda)*sin(2*pi*y(j)/rlambda)
         vx(i,j)=0.0
         vy(i,j)=0.0
         bx(i,j)=0.0
         by(i,j)=sqrt(4*pi*betai)
      enddo
      enddo

      return
      end
