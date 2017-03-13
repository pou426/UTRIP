c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),vx(ix,jx,kx)
      dimension vy(ix,jx,kx),vz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.4d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/real(ix-margin*2)
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

c-----------------------------------------------------------------------
c      dy,y

      dy0=1.d0/real(jx-margin*2)
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

c-----------------------------------------------------------------------
c      dz,z

      dz0=1.d0/real(kx-margin*2)
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
      ro1=0.125
      pr1=0.1

      pi = acos(-1.0d0)
      thini=60./180.*pi
      phini=30./180.*pi

      unitx=sin(thini)*cos(phini)
      unity=sin(thini)*sin(phini)
      unitz=cos(thini)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         if (unitx*x(i)+unity*y(j)+unitz*z(k).le.0.0d0) then
           ro(i,j,k)  = 1.
           pr(i,j,k)  = 1.
           vx(i,j,k)  = 0.
           vy(i,j,k)  = 0.
           vz(i,j,k)  = 0.
         else
           ro(i,j,k)  = ro1
           pr(i,j,k)  = pr1
           vx(i,j,k)  = 0.
           vy(i,j,k)  = 0.
           vz(i,j,k)  = 0.
         endif
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'phini',phini)


      
      return
      end
