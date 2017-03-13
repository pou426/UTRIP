c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,gx,gxm,gy,gym
     &            ,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx),gy(ix,jx),gym(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5./3.

      pi = acos(-1.0d0)

      xmax=2.d0
      xmin=0.4d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=(xmax-xmin)/real(ix-margin*2+1)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=xmin
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
c   gravitation
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         gx(i,j)=-1./x(i)**2
         gxm(i,j)=-1./(x(i)+dxm(i)/2)**2
         gy(i,j)=0.d0
         gym(i,j)=0.d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      amp=1.12d0
      amp=1.d0
      dr=0.05
      pr0=0.01

      do j=1,jx
      do i=1,ix
         ro(i,j)  = 1.d0+(amp-1.d0)*exp(- ((x(i)-1.d0)/dr)**2/2)
         pr(i,j)  = pr0/gm*ro(i,j)**gm
         vx(i,j) = 0.0
         vy(i,j) = 1.d0/sqrt(x(i))
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'amp',amp)
      call dacputparamd(mf_params,'dr',dr)


      
      return
      end
