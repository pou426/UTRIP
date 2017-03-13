c======================================================================|
      subroutine dx2dxi(dx,dxi,ix)
c======================================================================|

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxi(ix)

      do i=1,ix
         dxi(i) = 1.d0/dx(i)
      enddo

      return
      end
