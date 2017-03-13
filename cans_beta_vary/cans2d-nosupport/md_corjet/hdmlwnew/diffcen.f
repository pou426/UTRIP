c======================================================================|
      subroutine diffcen(un,u,dt,qxm,dx,dxm,qym,dy,dym,ix,jx)
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
      dimension un(ix,jx),u(ix,jx),qxm(ix,jx),qym(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)
c----------------------------------------------------------------------|
      call dx2dxi(dx,dxi,ix)
      call dx2dxi(dxm,dxim,ix)
      call dx2dxi(dy,dyi,jx)
      call dx2dxi(dym,dyim,jx)

      do j=1,jx
      do i=1,ix-1
        fx(i,j)=qxm(i,j)*dxim(i)*(u(i+1,j)-u(i,j))
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
        fy(i,j)=qym(i,j)*dyim(j)*(u(i,j+1)-u(i,j))
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
         un(i,j)=un(i,j)
     &        +dt*dxi(i)*(fx(i,j)-fx(i-1,j))
     &        +dt*dyi(j)*(fy(i,j)-fy(i,j-1))
      enddo
      enddo

      return
      end
