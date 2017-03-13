c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,omega0,qt0,dvy,yrgn,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

      omega0=1.0d-3
      qt0=1.5d0
      xrgn=2.d0
      yrgn=0.1d0
      zrgn=1.d0
      dvy=-qt0*omega0*xrgn

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xrgn/dble(ix-margin*2)
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

      dy0=yrgn/dble(jx-margin*2)
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

      dz0=zrgn/dble(kx-margin*2)
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
      betai=1.0d0/(1.d3)
      pr0=1.d-5
      b0=sqrt(8*pi*betai*pr0)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k) = 1.d0
         pr(i,j,k) = pr0
         bx(i,j,k) = 0.0d0
         by(i,j,k) = 0.0d0
         bz(i,j,k) = b0*(tanh((x(i)+0.2)/0.05)+1)/2.
     &               *(tanh((-x(i)+0.2)/0.05)+1)/2.
         vx(i,j,k) = 0.0d0
         vy(i,j,k) = -qt0*omega0*x(i)
         vz(i,j,k) = 0.0d0
      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      amp=1.d-2

      mdum=-1
      do k=1,kx
      do j=1,jx
      do i=1,ix
         rannum=2.d0*rangen(mdum)-1.d0
         pr(i,j,k) = pr(i,j,k)*(1+amp*rannum)
      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'omega0',omega0)
      call dacputparamd(mf_params,'qt0',qt0)
      call dacputparamd(mf_params,'xrgn',xrgn)
      call dacputparamd(mf_params,'yrgn',yrgn)
      call dacputparamd(mf_params,'zrgn',zrgn)



      return
      end
