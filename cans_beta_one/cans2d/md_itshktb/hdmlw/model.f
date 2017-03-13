c======================================================================|
      subroutine model(ro,vx,vy,cs2,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      cs2=1.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1./real(ix-margin*2-1)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=ix/2+1
      x(izero)=0.
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo


      dy0=1./real(jx-margin*2-1)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=jx/2+1
      y(jzero)=0.
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      ro0=1.
      ro1=0.125
      wtr=0.02

      pi = acos(-1.0d0)
      thini=60./180.*pi

      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'cs2',cs2)
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'ro1',ro1)

      
      return
      end
