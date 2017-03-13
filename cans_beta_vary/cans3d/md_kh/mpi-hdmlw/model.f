c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,margin,x,ix,y,jx,z,kx
     &           ,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      pi =      acos(-1.0d0)
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
      do kg=1,kgx
      do jg=1,jgx
      do ig=1,igx
        vrand=rangen(mdum)
        k=kg-kpe*(kx-2*margin)
        j=jg-jpe*(jx-2*margin)
        i=ig-ipe*(ix-2*margin)
        if (k.ge.1.and.k.le.kx .and.
     &      j.ge.1.and.j.le.jx .and.
     &      i.ge.1.and.i.le.ix) then
          rannum=2.d0*vrand-1.d0
          vz(i,j,k) = amp*rannum
        endif
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
