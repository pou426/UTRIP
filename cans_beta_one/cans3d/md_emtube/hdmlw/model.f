c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &   ,gx,gxm,gy,gym,gz,gzm,margin,x,ix,y,jx,z,kx
     &   ,rtube,ytube,ztube,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension gx(ix,jx,kx),gxm(ix,jx,kx)
      dimension gy(ix,jx,kx),gym(ix,jx,kx)
      dimension gz(ix,jx,kx),gzm(ix,jx,kx)

      dimension tem0(kx),pre0(kx),den0(kx),gzm0(kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.0d0/3.0d0
      zmax=25.0d0
      zmin=-22.5d0
      xmax=64.d0
      ymax=64.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=xmax/dble(ix-margin*2)
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

      dy0=ymax/dble(jx-margin*2)
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

      dz0=(zmax-zmin)/dble(kx-margin*2)
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=int(-zmin/dz0)+margin+1
      z(kzero)=dzm(kzero)/2.d0
      do k=kzero+1,kx
        z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
        z(k) = z(k+1)-dzm(k)
      enddo

c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      g0= -1.d0/gm

      do k=1,kx
        gzm0(k)=g0
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
         gx(i,j,k) =0.
         gxm(i,j,k)=0.
         gy(i,j,k) =0.
         gym(i,j,k)=0.
         gz(i,j,k) =g0
         gzm(i,j,k)=g0
      enddo
      enddo
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

c   temperature

      wtr=0.5d0
      tpho=1.0d0
      tcor=25.0d0
      dtem0=1.0d0

      zpho= 0.0d0
      ztr = 10.0d0
      zcor= 14.0d0

      dtem = dtem0*(gm-1.)/gm
      do k=1,kx
        if (z(k).lt.zpho) then
            tem0(k)=tpho-dtem*(z(k)-zpho)
        elseif (z(k).ge.zpho.and.z(k).lt.ztr) then
            tem0(k)=tpho
        elseif (z(k).ge.ztr.and.z(k).lt.zcor) then
            tem0(k)=tpho+(tcor-tpho)*(z(k)-ztr)/(zcor-ztr)
        elseif (z(k).ge.zcor) then
            tem0(k)=tcor
        endif
      enddo

c   density, pressure

      den0(kzero) = 1.d0
      pre0(kzero) = 1.d0/gm * tem0(kzero)
      do k=kzero+1,kx
         den0(k) = den0(k-1)
     &            *(tem0(k-1)+0.5*gm*gzm0(k-1)*dzm(k-1))
     &            /(tem0(k)  -0.5*gm*gzm0(k-1)*dzm(k-1))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo
      do k=kzero-1,1,-1
         den0(k) = den0(k+1)
     &             *(tem0(k+1)-0.5*gm*gzm0(k)*dzm(k))
     &             /(tem0(k)  +0.5*gm*gzm0(k)*dzm(k))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo

c   magnetic field

      rtube= 4.0d0
      ytube= 0.0d0
      ztube=-14.0d0
      btube=20.d0
      qtube=0.2d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         bx(i,j,k)= 0.d0
         by(i,j,k)= 0.d0
         bz(i,j,k)= 0.d0
         rr=sqrt((y(j)-ytube)**2+(z(k)-ztube)**2)
         if (rr.le.rtube) then
           bx(i,j,k)= btube/(1.0+(qtube*rr)**2.0)
           by(i,j,k)=-btube*qtube/(1.0+(qtube*rr)**2.0)*(z(k)-ztube)
           bz(i,j,k)= btube*qtube/(1.0+(qtube*rr)**2.0)*(y(j)-ytube)
         endif
      enddo
      enddo
      enddo

c  set all of aboves

      bxedge= btube/(1.0+(qtube*rtube)**2.0)
      btedge=-btube*qtube/(1.0+(qtube*rtube)**2.0)*rtube
      pmedge= (bxedge**2+btedge**2)/8.d0/pi

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k)=den0(k)
         pr(i,j,k)=pre0(k)
         rr=sqrt((y(j)-ytube)**2+(z(k)-ztube)**2)
         if (rr.le.rtube) then
           pr(i,j,k)=pr(i,j,k)-pmedge
         endif
         vx(i,j,k)=0.0d0
         vy(i,j,k)=0.0d0
         vz(i,j,k)=0.0d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ztr',ztr)


      return
      end
