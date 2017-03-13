c======================================================================|
      subroutine model(ro,pr,gm,rkap0
     &       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm
     &   ,margin,ix,jx,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dx(ix),xm(ix),dxm(ix)
      dimension y(jx),dy(jx),ym(jx),dym(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      rkap0=1.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=2.d0/dble(max(ix-margin*2,1))
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=ix/2+1
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

      dy0=2.d0/dble(max(jx-margin*2,1))
      do j=1,jx
         dym(j)=dy0
      enddo
      
      jzero=jx/2+1
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

      dz0=2.d0/dble(max(kx-margin*2,1))
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=kx/2+1
      if (kzero.lt.kx) then
        z(kzero)=dzm(kzero)/2.d0
      else
        kzero=1
        z(kzero)=0.d0
      endif
      kzero=2
      z(kzero)=0.d0
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo

      call grdrdy2(dz,zm,dzm,z,kx)

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
