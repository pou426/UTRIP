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
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1./real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=dx0/2
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=1./real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=jx/2+1
      y(jzero)=dy0/2
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      xshock=0.08d0
      xcont0=0.15d0
      wptb=0.024
      wtr=0.01


      do j=1,jx
      do i=1,ix
         xcont=xcont0-wptb*cos(4*pi*y(j))
         ro(i,j) = 1.23d0 
     &          + (0.17d0-1.23d0)*(1.d0+tanh((x(i)-xcont)/wtr))/2
         pr(i,j) = 1.01d0
         if (x(i).le.xshock) then
           ro(i,j) = 1.73d0
           pr(i,j) = 1.65d0
         endif
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)

      
      return
      end
