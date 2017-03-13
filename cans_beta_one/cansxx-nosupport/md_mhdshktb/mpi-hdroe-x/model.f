c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm,mu
     &       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm
     &       ,margin,ix,jx,kx,mf_params
     &       ,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dx(ix),xm(ix),dxm(ix)
      dimension y(jx),dy(jx),ym(jx),dym(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

      dimension xg(igx),dxg(igx),xmg(igx),dxmg(igx)
      dimension yg(jgx),dyg(jgx),ymg(jgx),dymg(jgx)
      dimension zg(kgx),dzg(kgx),zmg(kgx),dzmg(kgx)
      double precision mu

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      mu =4.d0*pi
      gm=2.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/dble(max(igx-margin*2,igx))
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

      call grdrdy2(dxg,xmg,dxmg,xg,igx)

      do i=1,ix
         ig=ipe*(ix-2*margin)+i
         x(i)=xg(ig)
         dx(i)=dxg(ig)
         xm(i)=xmg(ig)
         dxm(i)=dxmg(ig)
      enddo

c-----------------------------------------------------------------------
c      dy,y

      dy0=1.d0/dble(max(jgx-margin*2,jgx))
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

      call grdrdy2(dyg,ymg,dymg,yg,jgx)

      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         y(j)=yg(jg)
         dy(j)=dyg(jg)
         ym(j)=ymg(jg)
         dym(j)=dymg(jg)
      enddo

c-----------------------------------------------------------------------
c      dz,z

      dz0=1.d0/dble(max(kgx-margin*2,kgx))
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

      call grdrdy2(dzg,zmg,dzmg,zg,kgx)

      do k=1,kx
         kg=kpe*(kx-2*margin)+k
         z(k)=zg(kg)
         dz(k)=dzg(kg)
         zm(k)=zmg(kg)
         dzm(k)=dzmg(kg)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      b0=sqrt(2.d0*mu/gm)

c-- direction of background B-field

      thini=90.d0/180.d0*pi
      phini=0.d0/180.d0*pi

c     unitx=sin(thini)*cos(phini)
c     unity=sin(thini)*sin(phini)
c     unitz=cos(thini)

      unitx=1.d0
      unity=0.d0
      unitz=0.d0

c-- direction of kinked B-field

c     ebx0=1.d0
c     eby0=2.d0
c     ebz0=-(ebx0*unitx+eby0*unity)/unitz
c     ebb0=sqrt(ebx0**2+eby0**2+ebz0**2)
c     ebx=ebx0/ebb0
c     eby=eby0/ebb0
c     ebz=ebz0/ebb0
      ebx=0.d0
      eby=1.d0
      ebz=0.d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         if (unitx*x(i)+unity*y(j)+unitz*z(k).le.0.0d0) then
           ro(i,j,k)  = 1.d0
           pr(i,j,k)  = 1.d0
           bx(i,j,k)  = b0*0.75d0*unitx +b0*ebx
           by(i,j,k)  = b0*0.75d0*unity +b0*eby
           bz(i,j,k)  = b0*0.75d0*unitz +b0*ebz
         else
           ro(i,j,k)  = 0.125d0
           pr(i,j,k)  = 0.1d0
           bx(i,j,k)  = b0*0.75d0*unitx -b0*ebx
           by(i,j,k)  = b0*0.75d0*unity -b0*eby
           bz(i,j,k)  = b0*0.75d0*unitz -b0*ebz
         endif
         vx(i,j,k) = 0.d0
         vy(i,j,k) = 0.d0
         vz(i,j,k) = 0.d0
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
