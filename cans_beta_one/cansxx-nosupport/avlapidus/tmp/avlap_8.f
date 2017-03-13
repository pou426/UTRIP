c======================================================================|
      subroutine avlap_8
     &    (dan1,dan2,dan3,dan4,dan5,dan6,dan7,dan8
     &    ,da1,da2,da3,da4,da5,da6,da7,da8
     &    ,vx,vy,vz,dt,qav,vvmin,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)

      dimension dan1(ix,jx,kx),dan2(ix,jx,kx)
      dimension dan3(ix,jx,kx),dan4(ix,jx,kx)
      dimension dan5(ix,jx,kx),dan6(ix,jx,kx)
      dimension dan7(ix,jx,kx),dan8(ix,jx,kx)
      dimension da1(ix,jx,kx),da2(ix,jx,kx),da3(ix,jx,kx),da4(ix,jx,kx)
      dimension da5(ix,jx,kx),da6(ix,jx,kx),da7(ix,jx,kx),da8(ix,jx,kx)

      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension qvf(ix,jx,kx)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        qvf(i,j,k)=qav*dxm(i)*max(0.d0,abs(vx(i+1,j,k)-vx(i,j,k))-vvmin)
      enddo
      enddo
      enddo
      mdir=1
      call diffcen(dan1,da1,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan2,da2,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan3,da3,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan4,da4,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan5,da5,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan6,da6,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan7,da7,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan8,da8,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      endif

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        qvf(i,j,k)=qav*dym(j)*max(0.d0,abs(vy(i,j+1,k)-vy(i,j,k))-vvmin)
      enddo
      enddo
      enddo
      mdir=2
      call diffcen(dan1,da1,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan2,da2,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan3,da3,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan4,da4,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan5,da5,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan6,da6,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan7,da7,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan8,da8,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      endif

      if (mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        qvf(i,j,k)=qav*dzm(k)*max(0.d0,abs(vz(i,j,k+1)-vz(i,j,k))-vvmin)
      enddo
      enddo
      enddo
      mdir=3
      call diffcen(dan1,da1,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan2,da2,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan3,da3,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan4,da4,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan5,da5,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan6,da6,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan7,da7,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      call diffcen(dan8,da8,qvf
     &     ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
      endif

c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|

      return
      end
