c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,margin,x,ix,y,jx,z,kx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      pi =      acos(-1.0d0)
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
      ro0=1.0d0
      ro1=0.8d0
      u0=0.d0
      u1=0.5d0
      pr0=1.d0/gm
      pr1=1.d0/gm
      wtr=0.01d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k)=ro0+(ro1-ro0)*(1+tanh(z(k)/wtr))/2
         pr(i,j,k)=pr0+(pr1-pr0)*(1+tanh(z(k)/wtr))/2
         vx(i,j,k)=u0+(u1-u0)*(1+tanh(z(k)/wtr))/2
         vy(i,j,k)=0.d0
         vz(i,j,k)=0.d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      rlambda=1.d0/2.d0
      amp=1.d-2

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

      alpha0=ro0/(ro0+ro1)
      alpha1=ro1/(ro0+ro1)

      wkx=2.d0*pi/rlambda
      omegar=wkx*(alpha0*u0+alpha1*u1)
      omegai=sqrt(wkx*wkx*alpha0*alpha1*(u0-u1)*(u0-u1))

      do k=1,kx
      do j=1,jx
      do i=1,ix
       if (z(k).lt.0) then
         vre=amp*exp(wkx*z(k))*((wkx*u0-omegar)*cos(wkx*x(i))+
     &        omegai*sin(wkx*x(i)))
         vim=amp*exp(wkx*z(k))*((wkx*u0-omegar)*sin(wkx*x(i))-
     &        omegai*cos(wkx*x(i)))
         vz(i,j,k)=vz(i,j,k)+vre
         vx(i,j,k)=vx(i,j,k)-vim
         pr(i,j,k)=pr(i,j,k)
     &         -ro(i,j,k)/wkx*(omegai*vre-(wkx*u0-omegar)*vim)
       else
         vre=amp*exp(-wkx*z(k))*((wkx*u1-omegar)*cos(wkx*x(i))
     &        +omegai*sin(wkx*x(i)))
         vim=amp*exp(-wkx*z(k))*((wkx*u1-omegar)*sin(wkx*x(i))
     &        -omegai*cos(wkx*x(i)))
         vz(i,j,k)=vz(i,j,k)+vre
         vx(i,j,k)=vx(i,j,k)+vim
         pr(i,j,k)=pr(i,j,k)
     &        +ro(i,j,k)/wkx*(omegai*vre-(wkx*u1-omegar)*vim)
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
      call dacputparamd(mf_params,'pr1',pr1)
      call dacputparamd(mf_params,'vx0',vx0)
      call dacputparamd(mf_params,'vx1',vx1)


      
      return
      end
