c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,gx,gxm,gy,gym,gz,gzm
     &     ,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension gx(ix,jx,kx),gy(ix,jx,kx),gz(ix,jx,kx)
      dimension gxm(ix,jx,kx),gym(ix,jx,kx),gzm(ix,jx,kx)
      dimension gzmz(kx),roz(kx),prz(kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/(ix-margin*2)
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

      
      dy0=1.d0/(jx-margin*2)
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


      dz0=1.d0/(kx-margin*2)
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
c     gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do k=1,kx
        gzmz(k)=g0
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
         gx(i,j,k)=0.0
         gxm(i,j,k)=0.0
         gy(i,j,k)=0.0
         gym(i,j,k)=0.0
         gz(i,j,k)=g0
         gzm(i,j,k)=g0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     initial condition 
c----------------------------------------------------------------------|
      ro0=1.0d0
      ro1=4.00d0
      wtr=0.01d0

      pr1=ro1/gm
      prz(kx)=pr1
      do k=kx-1,1,-1
         if (z(k).le.0) then
           roz(k)=ro0
         else
           roz(k)=ro1
         endif
         roz(k)=ro0+(ro1-ro0)*(tanh(z(k)/wtr)+1)/2
         prz(k)=prz(k+1)-roz(k)*gzmz(k)*dzm(k)
      enddo

      do i=1,ix
      do j=1,jx
      do k=1,kx
         ro(i,j,k)=roz(k)
         pr(i,j,k)=prz(k)
         vx(i,j,k)=0.d0
         vy(i,j,k)=0.d0
         vz(i,j,k)=0.d0
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      rlambda=1.d0/2.d0
      amp=0.01d0

      if (.false.) then
      mdum=-1
      do k=1,kx
      do j=1,jx
      do i=1,ix
         rannum=2.d0*rangen(mdum)-1.d0
         vz(i,j,k) = amp*rannum
      enddo
      enddo
      enddo

      else

      wkx=2.d0*pi/rlambda
      omegai=sqrt(g0*wkx*(ro0-ro1)/(ro0+ro1))

      do k=1,kx
      do j=1,jx
      do i=1,ix
         if (z(k).le.0) then
           vv0=amp*exp(wkx*z(k))*omegai
           vx(i,j,k)= vv0*sin(wkx*x(i))
           vz(i,j,k)= vv0*cos(wkx*x(i))
           ro(i,j,k)=ro(i,j,k)*(1.d0+omegai*vz(i,j,k)/g0)
         else
           vv0=amp*exp(-wkx*z(k))*omegai
           vx(i,j,k)=-vv0*sin(wkx*x(i))
           vz(i,j,k)= vv0*cos(wkx*x(i))
           ro(i,j,k)=ro(i,j,k)*(1.d0+omegai*vz(i,j,k)/g0)
         endif
      enddo
      enddo
      enddo

      endif

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ro1',ro1)

      
      return
      end
