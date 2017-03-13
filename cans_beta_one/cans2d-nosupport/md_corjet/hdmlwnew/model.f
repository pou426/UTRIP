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

      ztr=8.0d0
      ymax=40.d0
      ymin=-4.d0

      pi = acos(-1.0d0)
      xmax= 40.d0
      xmin=-40.d0
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

      jzero=int(-ymin/dy0)+margin+1
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

      do j=1,jx
        gym0(j)=g0
      enddo

      do j=1,jx
      do i=1,ix
         gx(i,j)=0.0d0
         gxm(i,j)=0.0d0
         gy(i,j)=g0
         gym(i,j)=g0
      enddo
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

c   temperature

      wtr=0.5d0
      tcor=25.0d0
      tadg=2.0d0
      zcnv= 0.0d0

      do j=1,jx
            tem0(j)=1.d0
     &        -tadg*(gm-1.)/gm*(y(j)-zcnv)*0.5*(1-tanh(y(j)-zcnv))
     &        +0.5*(tcor-1.)*(tanh((y(j)-ztr)/wtr)+1.)
      enddo

c   plasma beta

      zf1=-2.0d0
      zf2=-0.0d0

      rbetaf=1/4.d0
      do j=1,jx
          af1 = (  tanh((y(j)-zf1)/0.5) + 1.)/2
          af2 = ( -tanh((y(j)-zf2)/0.5) + 1.)/2
          rbeta0(j) = rbetaf*af1*af2
      enddo

c   density, pressure

      den0(jzero) = 1.d0
      pre0(jzero) = 1.d0/gm * tem0(jzero)
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

c   magnetic field

      bxadd=-0.1d0
      do j=1,jx
          bmx0(j)  = sqrt(8*pi*pre0(j)*rbeta0(j))+bxadd
      enddo

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
      call dacputparamd(mf_params,'ztr',ztr)
      call dacputparamd(mf_params,'rbetaf',rbetaf)



      return
      end
