c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &   ,gx,gxm,gy,gym,gz,gzm,margin,x,ix,y,jx,z,kx
     &   ,rtube,ytube,ztube,mf_params
     &   ,igx,jgx,kgx,ipe,jpe,kpe)
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

      dimension tem0(kgx),pre0(kgx),den0(kgx),gzm0(kgx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

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

c-----------------------------------------------------------------------
c      dx,x

      dx0=xmax/dble(igx-margin*2)
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

      dy0=ymax/dble(jgx-margin*2)
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

      dz0=(zmax-zmin)/dble(kgx-margin*2)
      do k=1,kgx
         dzmg(k)=dz0
      enddo

      kzero=int(-zmin/dz0)+margin+1
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
c   gravitation
c----------------------------------------------------------------------|
      g0= -1.d0/gm

      do kg=1,kgx
        gzm0(kg)=g0
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
      do kg=1,kgx
        if (zg(kg).lt.zpho) then
            tem0(kg)=tpho-dtem*(zg(kg)-zpho)
        elseif (zg(kg).ge.zpho.and.zg(kg).lt.ztr) then
            tem0(kg)=tpho
        elseif (zg(kg).ge.ztr.and.zg(kg).lt.zcor) then
            tem0(kg)=tpho+(tcor-tpho)*(zg(kg)-ztr)/(zcor-ztr)
        elseif (zg(kg).ge.zcor) then
            tem0(kg)=tcor
        endif
      enddo

c   density, pressure

      den0(kzero) = 1.d0
      pre0(kzero) = 1.d0/gm * tem0(kzero)
      do k=kzero+1,kgx
         den0(k) = den0(k-1)
     &         *(tem0(k-1)+0.5*gm*gzm0(k-1)*dzmg(k-1))
     &         /(tem0(k)  -0.5*gm*gzm0(k-1)*dzmg(k-1))
         pre0(k) = pre0(kzero)
     &        *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
      enddo
      do k=kzero-1,1,-1
         den0(k) = den0(k+1)
     &             *(tem0(k+1)-0.5*gm*gzm0(k)*dzmg(k))
     &             /(tem0(k)  +0.5*gm*gzm0(k)*dzmg(k))
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
         kg=kpe*(kx-2*margin)+k
         ro(i,j,k)=den0(kg)
         pr(i,j,k)=pre0(kg)
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
