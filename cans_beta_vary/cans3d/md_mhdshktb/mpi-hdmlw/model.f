c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &       ,margin,x,ix,y,jx,z,kx,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=2.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/real(igx-margin*2)
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

      dy0=1.d0/real(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo

      jzero=jgx/2+1
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

      dz0=1.d0/real(kgx-margin*2)
      do k=1,kgx
         dzmg(k)=dz0
      enddo

      kzero=kgx/2+1
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
