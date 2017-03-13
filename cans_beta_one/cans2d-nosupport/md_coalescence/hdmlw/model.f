c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm
     &           ,et,etm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension et(ix,jx),etm(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      pi = acos(-1.0d0)
      xmax=2.d0*pi
      xmin=-xmax
      ymax=4.d0
      ymin=-ymin
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=(xmax-xmin)/dble(ix-margin*2)
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

      dy0=(ymax-ymin)/dble(jx-margin*2)
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
c   resistivity
c----------------------------------------------------------------------|
      etai=1.d-3
      do j=1,jx
      do i=1,ix
           et(i,j)=etai
           etm(i,j)=etai
      enddo
      enddo
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      betai=1.
      bb0=sqrt(8*pi*1/gm*betai)
      alpha=0.2

      do j=1,jx
      do i=1,ix
         denm=cosh(y(j))+alpha*cos(x(i))
         bx(i,j)=-bb0*sinh(y(j))/denm
         by(i,j)=-bb0*alpha*sin(x(i))/denm
         pr(i,j)=bb0**2/pi/8*(1-alpha**2)*(denm)**(-2)
         ro(i,j)=1.
         vx(i,j)=0.0
         vy(i,j)=0.0
      enddo
      enddo

      amp=1.d-2
      mdum=-1
      do j=1,jx
      do i=1,ix
         rannum=2.d0*rangen(mdum)-1.d0
         pr(i,j) = pr(i,j)*(1.d0+amp*rannum)
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'alpha',alpha)
      call dacputparamd(mf_params,'etai',etai)


      return
      end
