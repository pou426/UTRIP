c======================================================================|
      subroutine dx2ux(dx,dxm,ux0,ux1,ix)
c======================================================================|

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),ux0(ix),ux1(ix)

      do i=2,ix-1
         ux1(i)  = 0.5d0*dxm(i-1)/dx(i)
         ux0(i)  = 0.5d0*dxm(i)/dx(i)
      enddo

      return
      end
