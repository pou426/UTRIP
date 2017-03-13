c======================================================================|
      subroutine mlwartv(u,du,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
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
      dimension  dxi(ix),dxim(ix)
      dimension  dyi(jx),dyim(jx)
      dimension u(ix,jx),du(ix,jx),qx(ix,jx),qy(ix,jx)
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         du(i,j)=du(i,j)
     &        +dt*(dxi(i)*(qx(i,j)*dxim(i)*(u(i+1,j)-u(i,j))
     &                - qx(i-1,j)*dxim(i-1)*(u(i,j)-u(i-1,j))))
     &        +dt*(dyi(j)*(qy(i,j)*dyim(j)*(u(i,j+1)-u(i,j))
     &                - qy(i,j-1)*dyim(j-1)*(u(i,j)-u(i,j-1))))
      enddo
      enddo

      return
      end
