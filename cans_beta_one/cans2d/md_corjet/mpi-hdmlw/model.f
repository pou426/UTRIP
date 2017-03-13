c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm
     &           ,gx,gxm,gy,gym,margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx)
      dimension gy(ix,jx),gym(ix,jx)

      dimension tem0(jgx),pre0(jgx),den0(jgx),gym0(jgx)
      dimension bmx0(jgx),rbeta0(jgx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)

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

c-----------------------------------------------------------------------
c      dx,x

      dx0=(xmax-xmin)/dble(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=igx/2+1
      xg(izero)=dxmg(izero)/2.d0
      do i=izero+1,igx
         xg(i) = xg(i-1)+dxmg(i-1)
      enddo
      do i=izero-1,1,-1
         xg(i) = xg(i+1)-dxmg(i)
      enddo

      do i=1,ix
         ig=ipe*(ix-2*margin)+i
         x(i)=xg(ig)
         dxm(i)=dxmg(ig)
      enddo

c-----------------------------------------------------------------------
c      dy,y

      dy0=(ymax-ymin)/dble(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo
       
      jzero=int(-ymin/dy0)+margin+1
      yg(jzero)=dymg(jzero)/2.d0
      do j=jzero+1,jgx
         yg(j) = yg(j-1)+dymg(j-1)
      enddo
      do j=jzero-1,1,-1
         yg(j) = yg(j+1)-dymg(j)
      enddo

      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         y(j)=yg(jg)
         dym(j)=dymg(jg)
      enddo
c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do jg=1,jgx
        gym0(jg)=g0
      enddo

      do j=1,jx
      do i=1,ix
         gx(i,j)=0.0
         gxm(i,j)=0.0
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

      do jg=1,jgx
            tem0(jg)=1.d0
     &        -tadg*(gm-1.)/gm*(yg(jg)-zcnv)*0.5*(1-tanh(yg(jg)-zcnv))
     &        +0.5*(tcor-1.)*(tanh((yg(jg)-ztr)/wtr)+1.)
      enddo

c   plasma beta

      zf1=-2.0d0
      zf2=-0.0d0

      rbetaf=1/4.d0
      do jg=1,jgx
          af1 = (  tanh((yg(jg)-zf1)/0.5) + 1.)/2
          af2 = ( -tanh((yg(jg)-zf2)/0.5) + 1.)/2
          rbeta0(jg) = rbetaf*af1*af2
      enddo

c   density, pressure

      den0(jzero) = 1.d0
      pre0(jzero) = 1.d0/gm * tem0(jzero)
      do jg=jzero+1,jgx
         den0(jg) = den0(jg-1)
     &    *((1+rbeta0(jg-1))*tem0(jg-1)+0.5*gm*gym0(jg-1)*dymg(jg-1))
     &    /((1+rbeta0(jg)  )*tem0(jg)  -0.5*gm*gym0(jg-1)*dymg(jg-1))
         pre0(jg) = pre0(jzero)
     &        *(den0(jg)/den0(jzero))*(tem0(jg)/ tem0(jzero))
      enddo
      do jg=jzero-1,1,-1
         den0(jg) = den0(jg+1)
     &    *((1+rbeta0(jg+1))*tem0(jg+1)-0.5*gm*gym0(jg)*dymg(jg))
     &    /((1+rbeta0(jg)  )*tem0(jg)  +0.5*gm*gym0(jg)*dymg(jg))
         pre0(jg) = pre0(jzero)
     &        *(den0(jg)/den0(jzero))*(tem0(jg)/ tem0(jzero))
      enddo

c   magnetic field

      bxadd=-0.1d0
      do jg=1,jgx
          bmx0(jg)  = sqrt(8*pi*pre0(jg)*rbeta0(jg))+bxadd
      enddo


c  set all of aboves 

      do j=1,jx
      do i=1,ix
         jg=jpe*(jx-2*margin)+j
         ro(i,j)=den0(jg)
         pr(i,j)=pre0(jg)
         vx(i,j)=0.0d0
         vy(i,j)=0.0d0
         bx(i,j)=bmx0(jg)
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
