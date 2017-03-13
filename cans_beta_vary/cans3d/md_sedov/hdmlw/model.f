c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi=acos(-1.0d0)
      gm=5.d0/3.d0

      xmin=0.04d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=xmin+dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=2.d0*pi/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=margin+1
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
       
      kzero=margin+1
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
      prism=1.d-8
      wexp=0.1d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k) = 1.0d0
         vx(i,j,k) = 0.0d0
         vy(i,j,k) = 0.0d0
         vz(i,j,k) = 0.0d0
         ss=sqrt(x(i)**2+z(k)**2)
         pr(i,j,k)  = prism+(1./gm-prism)*exp(-(ss/wexp)**2)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'prism',prism)
      call dacputparamd(mf_params,'wexp',wexp)


      
      return
      end
