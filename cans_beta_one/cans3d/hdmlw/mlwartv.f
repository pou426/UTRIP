c======================================================================|
      subroutine mlwartv(u,du,dt
     &             ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
c======================================================================|
c
c NAME  mlwartv
c
c PURPOSE
c    apply artificial viscosity
c
c INPUTS & OUTPUTS
c    du(ix,jx,kx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    u(ix,jx,kx) : [double] basic variables
c    qx(ix,jx,kx) : [double] coefficients of artificial viscosity
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension  dxi(ix),dxim(ix)
      dimension  dyi(jx),dyim(jx)
      dimension  dzi(kx),dzim(kx)
      dimension  u(ix,jx,kx),du(ix,jx,kx),qx(ix,jx,kx),qy(ix,jx,kx)
      dimension qz(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         du(i,j,k)=du(i,j,k)
     &        +dt*(dxi(i)*(qx(i,j,k)*dxim(i)*(u(i+1,j,k)-u(i,j,k))
     &                - qx(i-1,j,k)*dxim(i-1)*(u(i,j,k)-u(i-1,j,k))))
     &        +dt*(dyi(j)*(qy(i,j,k)*dyim(j)*(u(i,j+1,k)-u(i,j,k))
     &                - qy(i,j-1,k)*dyim(j-1)*(u(i,j,k)-u(i,j-1,k))))
     &        +dt*(dzi(k)*(qz(i,j,k)*dzim(k)*(u(i,j,k+1)-u(i,j,k))
     &                - qz(i,j,k-1)*dzim(k-1)*(u(i,j,k)-u(i,j,k-1))))
      enddo
      enddo
      enddo

      return
      end
