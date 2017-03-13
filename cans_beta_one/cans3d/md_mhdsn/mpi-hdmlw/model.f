c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm,margin
     &                ,x,ix,y,jx,z,kx
     &           ,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi=acos(-1.0d0)
      gm=5.d0/3.d0

      xmin=0.04d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/real(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=margin+1
      xg(izero)=xmin+dxmg(izero)/2.d0
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

      dy0=2.d0*pi/real(jgx-margin*2)
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

      dz0=1.d0/real(kgx-margin*2)
      do k=1,kgx
         dzmg(k)=dz0
      enddo

      kzero=margin+1
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
      prism=1.d-8
      wexp=0.1d0

      betai=1.0d6
      b0=sqrt(prism*8*pi*betai)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k) = 1.0d0
         vx(i,j,k) = 0.0d0
         vy(i,j,k) = 0.0d0
         vz(i,j,k) = 0.0d0
         ss=sqrt(x(i)**2+z(k)**2)
         pr(i,j,k) = prism+(1.d0/gm-prism)*exp(-(ss/wexp)**2)
         bx(i,j,k) = 0.0d0
         by(i,j,k) = 0.0d0
         bz(i,j,k) = b0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'wexp',wexp)
      call dacputparamd(mf_params,'prism',prism)
      call dacputparamd(mf_params,'betai',betai)
      call dacputparamd(mf_params,'gm',gm)

      
      return
      end
