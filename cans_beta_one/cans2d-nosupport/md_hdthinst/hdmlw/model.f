c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)
      parameter(kmax=2)
      dimension phase(-kmax:kmax,-kmax:kmax)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      xmax=0.3d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=xmax/real(ix-margin*2)
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

      dy0=xmax/real(jx-margin*2)
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

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      ro0=2.93289
      pr0=2.21862

      pi = acos(-1.0d0)

      iseed=10000
      do kx=-kmax,kmax
      do ky=-kmax,kmax
         phase(kx,ky) = rangen(iseed)
      enddo
      enddo


      amp=0.15d0

      do j=1,jx
      do i=1,ix

        cosi = 0.0d0
        sine = 0.0d0

        do kx = -kmax, kmax
        do ky = -kmax, kmax
          wkxky=(kx*x(i) + ky*y(j))/xmax
          wkxky=2*pi*( wkxky + phase(kx,ky))
          wkabs=(kx**2+ky**2)**.5
          if(wkabs.ne.0)then
             cosi = cosi + cos(wkxky)
             sine = sine + sin(wkxky)
          endif
        enddo
        enddo

        cosi = amp * cosi / real(2*kmax+1)
        sine = amp * sine / real(2*kmax+1)

         ro(i,j) = ro0* ( 1.0d0 + cosi )
         pr(i,j) = pr0
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|

      
      return
      end
