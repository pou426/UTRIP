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
      gm=5./3.

      pi = acos(-1.0d0)

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1./real(ix-margin*2+1)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      izero=1
      x(izero)=0.5d0*dx0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=2.d0*pi/real(jx-margin*2+1)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=margin+1
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
      prism=1.e-8
      wexp=0.02

      do j=1,jx
      do i=1,ix
         ro(i,j)  = 1.
         vx(i,j) = 0.0
         vy(i,j) = 0.0
         ss=abs(x(i))
         pr(i,j)  = prism+(1.-prism)*exp(-(ss/wexp)**2)
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'prism',prism)
      call dacputparamd(mf_params,'wexp',wexp)

      
      return
      end
