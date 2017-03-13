c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=1.d0/dble(ix-margin*2)
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

      dy0=1.d0/dble(jx-margin*2)
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
      betai=0.5
      b0=sqrt(8*pi/gm*betai)

      thini=60./180.*pi
      wexp=0.02
      amp=0.1d0

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.
         ss=sqrt(x(i)**2+y(j)**2)
         pr(i,j) = 1/gm+amp/gm*exp(-(ss/wexp)**2)
         bx(i,j) = b0*cos(thini)
         by(i,j) = b0*sin(thini)
c        bx(i,j) = 0.0
c        by(i,j) = b0
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'thini',thini)


      return
      end
