c======================================================================|
      subroutine diffcen(un,u
     &    ,dt,qxm,dx,dxm,qym,dy,dym,qzm,dz,dzm,ix,jx,kx,mfdim)
c======================================================================|
c
c NAME  mlwartv
c
c PURPOSE
c    apply artificial viscosity
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    u(ix,jx) : [double] basic variables
c    qx(ix,jx) : [double] coefficients of artificial viscosity
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension  dx(ix),dxm(ix),dxi(ix),dxim(ix)
      dimension  dy(jx),dym(jx),dyi(jx),dyim(jx)
      dimension  dz(kx),dzm(kx),dzi(kx),dzim(kx)
      dimension un(ix,jx,kx),u(ix,jx,kx)
      dimension qxm(ix,jx,kx),qym(ix,jx,kx),qzm(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
c----------------------------------------------------------------------|

      if (mfdim(1).eq.1) then
      call dx2dxi(dx,dxi,ix)
      call dx2dxi(dxm,dxim,ix)
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        fx(i,j,k)=qxm(i,j,k)*dxim(i)*(u(i+1,j,k)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx
      do i=2,ix
         un(i,j,k)=un(i,j,k)
     &        +dt*dxi(i)*(fx(i,j,k)-fx(i-1,j,k))
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      call dx2dxi(dy,dyi,jx)
      call dx2dxi(dym,dyim,jx)
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        fy(i,j,k)=qym(i,j,k)*dyim(j)*(u(i,j+1,k)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx
      do i=1,ix
         un(i,j,k)=un(i,j,k)
     &        +dt*dyi(j)*(fy(i,j,k)-fy(i,j-1,k))
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      call dx2dxi(dz,dzi,kx)
      call dx2dxi(dzm,dzim,kx)
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        fz(i,j,k)=qzm(i,j,k)*dzim(k)*(u(i,j,k+1)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx
      do i=1,ix
         un(i,j,k)=un(i,j,k)
     &        +dt*dzi(k)*(fz(i,j,k)-fz(i,j,k-1))
      enddo
      enddo
      enddo
      endif

      return
      end
