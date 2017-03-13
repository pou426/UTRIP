c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,omega0,qt0,dvy,yrgn,margin,x,ix,y,jx,z,kx,mf_params
     &    ,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension dzm(kx),z(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

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
c-----------------------------------------------------------------------
c      dx,x

      dx0=xrgn/dble(igx-margin*2)
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

      dy0=yrgn/dble(jgx-margin*2)
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

      dz0=zrgn/dble(kgx-margin*2)
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
      do kg=1,kgx
      do jg=1,jgx
      do ig=1,igx
        rannum=2.d0*rangen(mdum)-1.d0
        k=kg-kpe*(kx-2*margin)
        j=jg-jpe*(jx-2*margin)
        i=ig-ipe*(ix-2*margin)
        if (k.ge.1.and.k.le.kx .and.
     &      j.ge.1.and.j.le.jx .and.
     &      i.ge.1.and.i.le.ix) then
          pr(i,j,k) = pr(i,j,k)*(1+amp*rannum)
        endif
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
