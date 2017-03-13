c======================================================================|
      subroutine getfby(fx,fy,ez,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ez(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      do j=1,jx
      do i=1,ix
         fx(i,j)= -ez(i,j)
         fy(i,j)= 0.d0
      enddo
      enddo

      return
      end
