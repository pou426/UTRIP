c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.4d0

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
      ro0=1.
      ro1=0.125
      pr0=1.
      pr1=0.1

      pi = acos(-1.0d0)
      thini=60./180.*pi

      wtr=0.02


      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)


      
      return
      end
