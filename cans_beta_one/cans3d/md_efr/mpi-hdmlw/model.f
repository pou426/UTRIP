c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &           ,gx,gxm,gy,gym,gz,gzm,margin,x,ix,y,jx,z,kx
     &           ,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension gx(ix,jx,kx),gxm(ix,jx,kx)
      dimension gy(ix,jx,kx),gym(ix,jx,kx)
      dimension gz(ix,jx,kx),gzm(ix,jx,kx)

      dimension tem0(kgx),pre0(kgx),den0(kgx),gzm0(kgx)
      dimension bmx0(kgx),rbeta0(kgx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      gm=5.d0/3.d0
      ztr=10.d0
      zmax=40.d0
      zmin=-4.d0

      xmax=40.d0
      ymax=2.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=xmax/real(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=margin+1
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

      dy0=ymax/real(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo
       
      jzero=margin+1
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

c-----------------------------------------------------------------------
c      dz,z

      dz0=(zmax-zmin)/real(kgx-margin*2)
      do k=1,kgx
         dzmg(k)=dz0
      enddo

      kzero=int(-zmin/dz0)+margin+1
      zg(kzero)=dzmg(kzero)/2.d0
      do k=kzero+1,kgx
         zg(k) = zg(k-1)+dzmg(k-1)
      enddo
      do k=kzero-1,1,-1
         zg(k) = zg(k+1)-dzmg(k)
      enddo

      do k=1,kx
         kg=kpe*(kx-2*margin)+k
         z(k)=zg(kg)
         dzm(k)=dzmg(kg)
      enddo
c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do kg=1,kgx
        gzm0(kg)=g0
      enddo


      do k=1,kx
      do j=1,jx
      do i=1,ix
         gx(i,j,k) =0.
         gxm(i,j,k)=0.
         gy(i,j,k) =0.
         gym(i,j,k)=0.
         gz(i,j,k) =g0
         gzm(i,j,k)=g0
      enddo
      enddo
      enddo
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

c   temperature

      wtr=0.5d0
      tcor=25.d0
      tadg=2.d0
      zcnv= 0.d0

      do k=1,kgx
            tem0(k)=1.
     &        -tadg*(gm-1.)/gm*(zg(k)-zcnv)*0.5*(1-tanh(zg(k)-zcnv))
     &        +0.5*(tcor-1.)*(tanh((zg(k)-ztr)/wtr)+1.)
      enddo

c   plasma beta

      zf1=-2.0d0
      zf2=-0.0d0

      rbetaf=1/4.d0
      do k=1,kgx
          af1 = (  tanh((zg(k)-zf1)/0.5) + 1.)/2
          af2 = ( -tanh((zg(k)-zf2)/0.5) + 1.)/2
          rbeta0(k) = rbetaf*af1*af2
      enddo

c   density, pressure

      den0(kzero) = 1.d0
      pre0(kzero) = 1.d0/gm * tem0(kzero)
      do k=kzero+1,kgx
         den0(k) = den0(k-1)
     &         *((1+rbeta0(k-1))*tem0(k-1)+0.5*gm*gzm0(k-1)*dzmg(k-1))
     &         /((1+rbeta0(k)  )*tem0(k)  -0.5*gm*gzm0(k-1)*dzmg(k-1))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo
      do k=kzero-1,1,-1
         den0(k) = den0(k+1)
     &             *((1+rbeta0(k+1))*tem0(k+1)-0.5*gm*gzm0(k)*dzmg(k))
     &             /((1+rbeta0(k)  )*tem0(k)  +0.5*gm*gzm0(k)*dzmg(k))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo

c   magnetic field

      do k=1,kgx
          bmx0(k)  = sqrt(8*pi*pre0(k)*rbeta0(k))
      enddo

c  set all of aboves

      do k=1,kx
      do j=1,jx
      do i=1,ix
         kg=kpe*(kx-2*margin)+k
         ro(i,j,k)=den0(kg)
         pr(i,j,k)=pre0(kg)
         vx(i,j,k)=0.0d0
         vy(i,j,k)=0.0d0
         vz(i,j,k)=0.0d0
         bx(i,j,k)=bmx0(kg)
         by(i,j,k)=0.0d0
         bz(i,j,k)=0.0d0
      enddo
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
