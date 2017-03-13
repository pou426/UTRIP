c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,gm,gx,gxm,gy,gym,gz,gzm
     &     ,margin,x,ix,y,jx,z,kx
     &           ,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
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
      dimension gzmz(kgx),roz(kgx),prz(kgx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

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
c     gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do k=1,kgx
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
      prz(kgx)=pr1
      do kg=kgx-1,1,-1
         roz(kg)=ro0+(ro1-ro0)*(tanh(zg(kg)/wtr)+1)/2
         prz(kg)=prz(kg+1)-roz(kg)*gzmz(kg)*dzmg(kg)
      enddo

      do i=1,ix
      do j=1,jx
      do k=1,kx
         kg=kpe*(kx-2*margin)+k
         ro(i,j,k)=roz(kg)
         pr(i,j,k)=prz(kg)
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

      wkx=2.0*pi/rlambda
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
