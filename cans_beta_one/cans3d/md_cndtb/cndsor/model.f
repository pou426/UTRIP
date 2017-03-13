c======================================================================|
      subroutine model(ro,pr,gm,rkap0,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      rkap0=1.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=2.d0/real(ix-margin*2)
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

      dy0=2.d0/real(jx-margin*2)
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

      dz0=2.d0/real(kx-margin*2)
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=kx/2+1
      z(kzero)=dzm(kzero)/2.d0
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo
c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      thini=60.d0/180.d0*pi
      phini=30.d0/180.d0*pi

      unitx=sin(thini)*cos(phini)
      unity=sin(thini)*sin(phini)
      unitz=cos(thini)

      wexp=0.3d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k) = 1.d0
         s=unitx*x(i)+unity*y(j)+unitz*z(k)
         pr(i,j,k) = 1/gm*exp(-(s/wexp)**2)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'phini',phini)

      
      return
      end
