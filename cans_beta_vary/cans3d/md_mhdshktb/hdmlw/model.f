c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &       ,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=2.d0

      pi = acos(-1.0d0)
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
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
      b0=sqrt(4*pi)

      thini=60./180.*pi
      phini=30./180.*pi

      unitx=sin(thini)*cos(phini)
      unity=sin(thini)*sin(phini)
      unitz=cos(thini)

      ebx0=1.
      eby0=2.
      ebz0=-(ebx0*unitx+eby0*unity)/unitz
      ebb0=sqrt(ebx0**2+eby0**2+ebz0**2)
      ebx=ebx0/ebb0
      eby=eby0/ebb0
      ebz=ebz0/ebb0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         if (unitx*x(i)+unity*y(j)+unitz*z(k).le.0.0d0) then
           ro(i,j,k)  = 1.
           pr(i,j,k)  = 1.
           bx(i,j,k)  = b0*0.75*unitx +b0*ebx
           by(i,j,k)  = b0*0.75*unity +b0*eby
           bz(i,j,k)  = b0*0.75*unitz +b0*ebz
         else
           ro(i,j,k)  = 0.125
           pr(i,j,k)  = 0.1
           bx(i,j,k)  = b0*0.75*unitx -b0*ebx
           by(i,j,k)  = b0*0.75*unity -b0*eby
           bz(i,j,k)  = b0*0.75*unitz -b0*ebz
         endif
         vx(i,j,k) = 0.0
         vy(i,j,k) = 0.0
         vz(i,j,k) = 0.0
      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'thini',thini)
      call dacputparamd(mf_params,'phini',phini)
      call dacputparamd(mf_params,'ebx',ebx)
      call dacputparamd(mf_params,'eby',eby)
      call dacputparamd(mf_params,'ebz',ebz)


      return
      end
