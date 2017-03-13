c======================================================================|
      subroutine avlap_6
     &    (dan1,dan2,dan3,dan4,dan5,dan6
     &    ,da1,da2,da3,da4,da5,da6
     &    ,vx,vy,dt,qav,vvmin,dx,dxm,dy,dym,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)

      dimension dan1(ix,jx),dan2(ix,jx),dan3(ix,jx),dan4(ix,jx)
      dimension dan5(ix,jx),dan6(ix,jx)
      dimension da1(ix,jx),da2(ix,jx),da3(ix,jx),da4(ix,jx)
      dimension da5(ix,jx),da6(ix,jx)

      dimension vx(ix,jx),vy(ix,jx)
      dimension qxm(ix,jx)
      dimension qym(ix,jx)

c-------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
      vvmin=1.d-4

      do j=1,jx-1
      do i=1,ix-1
         qxm(i,j)=qav*dxm(i)*max(0.d0,abs(vx(i+1,j)-vx(i,j))-vvmin)
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
         qym(i,j)=qav*dym(j)*max(0.d0,abs(vy(i,j+1)-vy(i,j))-vvmin)
      enddo
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call diffcen(dan1,da1,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
      call diffcen(dan2,da2,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
      call diffcen(dan3,da3,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
      call diffcen(dan4,da4,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
      call diffcen(dan5,da5,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
      call diffcen(dan6,da6,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)

      return
      end
