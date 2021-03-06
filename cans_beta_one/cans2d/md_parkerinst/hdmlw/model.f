c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm
     &           ,gx,gxm,gy,gym,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx)
      dimension gy(ix,jx),gym(ix,jx)

      dimension tem0(jx),pre0(jx),den0(jx),gym0(jx)
      dimension bmx0(jx),rbeta0(jx)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|

      gm=5.d0/3.d0

      pi = acos(-1.0d0)

      xmax=50.d0
      xmin=-xmax
      ymax=25.d0
      ymin=-ymax
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
c   gravitation
c----------------------------------------------------------------------|
       g0= -1.0d0/gm
       wg0=0.5

       do j=1,jx
        gym0(j)=g0 * tanh((y(j)+0.5*dym(j))/wg0)
       enddo

       do j=1,jx
       do i=1,ix
        gx(i,j)=0.0
        gxm(i,j)=0.0
        gy(i,j)=g0 * tanh(y(j)/wg0)
        gym(i,j)=g0 * tanh((y(j)+0.5*dym(j))/wg0)
       enddo
       enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

c-------------------------
c   temperature

       ytr=6.0
       wtr=0.5
       tcor=25.0

       do j=1,jx
            dycen=abs(y(j))
            tem0(j)=1.
     &        +0.5*(tcor-1.)*(tanh((dycen-ytr)/wtr)+1.)
       enddo


c-------------------------
c   plasma beta

      yf1=-2.5
      yf2= 2.5
      rbetaf=1.0/0.2
      wf0=0.5

      do j=1,jx
          af1 = (  tanh((y(j)-yf1)/wf0) + 1.)/2
          af2 = ( -tanh((y(j)-yf2)/wf0) + 1.)/2
          rbeta0(j) = abs(rbetaf*af1*af2)
      enddo

c-------------------------
c   density, pressure

      den0(jzero) = 1.
      pre0(jzero) = 1./gm * tem0(jzero)

      do j=jzero+1,jx
         den0(j) = den0(j-1)
     &            *((1+rbeta0(j-1))*tem0(j-1)+0.5*gm*gym0(j-1)*dym(j-1))
     &            /((1+rbeta0(j)  )*tem0(j)  -0.5*gm*gym0(j-1)*dym(j-1))
         pre0(j) = pre0(jzero)
     &        *(den0(j)/den0(jzero))*(tem0(j)/ tem0(jzero))
      enddo
      do j=jzero-1,1,-1
         den0(j) = den0(j+1)
     &             *((1+rbeta0(j+1))*tem0(j+1)-0.5*gm*gym0(j)*dym(j))
     &             /((1+rbeta0(j)  )*tem0(j)  +0.5*gm*gym0(j)*dym(j))
         pre0(j) = pre0(jzero)
     &        *(den0(j)/den0(jzero))*(tem0(j)/ tem0(jzero))
      enddo


c-------------------------
c   magnetic field

       do j=1,jx
          bmx0(j)  = sqrt(8*pi*pre0(j)*rbeta0(j))
       enddo

c-------------------------
c  set all of aboves 

      do j=1,jx
      do i=1,ix
         ro(i,j)=den0(j)
         pr(i,j)=pre0(j)
         vx(i,j)=0.0d0
         vy(i,j)=0.0d0
         bx(i,j)=bmx0(j)
         by(i,j)=0.0d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'rbetaf',rbetaf)

      return
      end

