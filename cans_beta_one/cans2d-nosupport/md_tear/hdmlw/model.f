c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &           ,et,etm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)
      dimension et(ix,jx),etm(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|

      gm=5./3.

      pi = acos(-1.0d0)
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      alpha=0.5
      xmax=2*pi/alpha
      dx0=xmax/real(ix-margin*2)
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

      ymax=6.d0
      dy0=ymax/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin+1
      y(jzero)=dy0/2
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      betai=5.d0/6.d0
      bb0=sqrt(betai*pi*8/gm)

      do j=1,jx
      do i=1,ix
         bx(i,j)=bb0*tanh(y(j))
         by(i,j)=0.0d0
         bz(i,j)=bb0/cosh(y(j))
         pr(i,j)=1/gm
         ro(i,j)=1.0d0
         vx(i,j)=0.0d0
         vy(i,j)=0.0d0
         vz(i,j)=0.0d0
      enddo
      enddo


c----------------------------------------------------------------------|
c   resistivity
c----------------------------------------------------------------------|
      rm=1.d2
      etai=2*betai/rm/gm

      do j=1,jx
      do i=1,ix
        et(i,j)=etai
        etm(i,j)=etai
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      amp=0.01

      do j=1,jx
      do i=1,ix
        vx(i,j) = vx(i,j) -amp*sin(2*pi*x(i)/xmax)/cosh(y(j))
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'betai',betai)


      return
      end
