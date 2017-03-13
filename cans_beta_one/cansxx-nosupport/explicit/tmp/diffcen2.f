c======================================================================|
      subroutine diffcen(un,u,qvf
     &    ,dt,dx,dxm,dy,dym,dz,dzm,ix,jx,kx,mdir)
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
      dimension  dx(ix),dxm(ix),dxi(ix),dxim(ix)
      dimension  dy(jx),dym(jx),dyi(jx),dyim(jx)
      dimension  dz(kx),dzm(kx),dzi(kx),dzim(kx)
      dimension un(ix,jx,kx),u(ix,jx,kx)
      dimension qvf(ix,jx,kx)
      dimension f(ix,jx,kx)
c----------------------------------------------------------------------|

      if (mdir.eq.1) then
      call dx2dxi(dx,dxi,ix)
      call dx2dxi(dxm,dxim,ix)
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        f(i,j,k)=-qvf(i,j,k)*dxim(i)*(u(i+1,j,k)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx
      do i=2,ix
         un(i,j,k)=un(i,j,k)
     &        -dt*dxi(i)*(f(i,j,k)-f(i-1,j,k))
      enddo
      enddo
      enddo
      endif

      if (mdir.eq.2) then
      call dx2dxi(dy,dyi,jx)
      call dx2dxi(dym,dyim,jx)
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        f(i,j,k)=-qvf(i,j,k)*dyim(j)*(u(i,j+1,k)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx
      do i=1,ix
         un(i,j,k)=un(i,j,k)
     &        -dt*dyi(j)*(f(i,j,k)-f(i,j-1,k))
      enddo
      enddo
      enddo
      endif

      if (mdir.eq.3) then
      call dx2dxi(dz,dzi,kx)
      call dx2dxi(dzm,dzim,kx)
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        f(i,j,k)=-qvf(i,j,k)*dzim(k)*(u(i,j,k+1)-u(i,j,k))
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx
      do i=1,ix
         un(i,j,k)=un(i,j,k)
     &        -dt*dzi(k)*(f(i,j,k)-f(i,j,k-1))
      enddo
      enddo
      enddo
      endif

      return
      end
