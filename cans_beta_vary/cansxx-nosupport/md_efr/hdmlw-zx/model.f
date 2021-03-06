c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm,mu,gz,gzm
     &       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm
     &       ,margin,ix,jx,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dx(ix),xm(ix),dxm(ix)
      dimension y(jx),dy(jx),ym(jx),dym(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension gz(ix,jx,kx),gzm(ix,jx,kx)
      double precision mu

      dimension tem0(kx),pre0(kx),den0(kx),gzm0(kx)
      dimension bmx0(kx),rbeta0(kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      pi = acos(-1.0d0)
      mu =4.d0*pi
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      xmax=40.d0
      ymax=2.d0
      ztr=10.d0
      zmax=40.d0
      zmin=-4.d0
c----------------------------------------------------------------------|
      dx0=xmax/dble(max(ix-margin*2,1))
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      if (izero.lt.ix) then
        x(izero)=dxm(izero)/2.d0
      else
        izero=1
        x(izero)=0.d0
      endif
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      call grdrdy2(dx,xm,dxm,x,ix)

      dy0=ymax/dble(max(jx-margin*2,1))
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin+1
      if (jzero.lt.jx) then
        y(jzero)=dym(jzero)/2.d0
      else
        jzero=1
        y(jzero)=0.d0
      endif
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

      call grdrdy2(dy,ym,dym,y,jx)

      dz0=(zmax-zmin)/dble(max(kx-margin*2,1))
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=int(-zmin/dz0)+margin+1
      if (kzero.lt.kx) then
        z(kzero)=0.d0
      else
        kzero=1
        z(kzero)=0.d0
      endif
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo

      call grdrdy2(dz,zm,dzm,z,kx)

c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do k=1,kx
        gzm0(k)=g0
      enddo


      do k=1,kx
      do j=1,jx
      do i=1,ix
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
      tcor=25.0d0
      tadg=2.0d0
      zcnv= 0.0d0

      do k=1,kx
            tem0(k)=1.
     &        -tadg*(gm-1.)/gm*(z(k)-zcnv)*0.5*(1-tanh(z(k)-zcnv))
     &        +0.5*(tcor-1.)*(tanh((z(k)-ztr)/wtr)+1.)
      enddo

c   plasma beta

      zf1=-2.0d0
      zf2=-0.0d0

      rbetaf=1/4.d0
      do k=1,kx
          af1 = (  tanh((z(k)-zf1)/0.5) + 1.)/2
          af2 = ( -tanh((z(k)-zf2)/0.5) + 1.)/2
          rbeta0(k) = rbetaf*af1*af2
      enddo

c   density, pressure

      den0(kzero) = 1.d0
      pre0(kzero) = 1.d0/gm * tem0(kzero)
      do k=kzero+1,kx
         den0(k) = den0(k-1)
     &         *((1+rbeta0(k-1))*tem0(k-1)+0.5*gm*gzm0(k-1)*dzm(k-1))
     &         /((1+rbeta0(k)  )*tem0(k)  -0.5*gm*gzm0(k-1)*dzm(k-1))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo
      do k=kzero-1,1,-1
         den0(k) = den0(k+1)
     &             *((1+rbeta0(k+1))*tem0(k+1)-0.5*gm*gzm0(k)*dzm(k))
     &             /((1+rbeta0(k)  )*tem0(k)  +0.5*gm*gzm0(k)*dzm(k))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo

c   magnetic field

      do k=1,kx
          bmx0(k)  = sqrt(8*pi*pre0(k)*rbeta0(k))
      enddo

c  set all of aboves

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k)=den0(k)
         pr(i,j,k)=pre0(k)
         vx(i,j,k)=0.0d0
         vy(i,j,k)=0.0d0
         vz(i,j,k)=0.0d0
         bx(i,j,k)=bmx0(k)
         by(i,j,k)=0.0d0
         bz(i,j,k)=0.0d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      amp=0.05d0
      xptb=20.d0
      zptb1=-2.d0
      zptb2= 0.d0
      wptb=0.5d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
        vz(i,j,k) = vz(i,j,k)
     &    +amp*cos(2.d0*pi*x(i)/xptb)
     &    *0.5d0*(tanh((x(i)+0.75d0*xptb)/wptb)
     &           -tanh((x(i)-0.75d0*xptb)/wptb))
     &    *0.5d0*(tanh((z(k)-zptb1)/wptb)
     &           -tanh((z(k)-zptb2)/wptb))
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
