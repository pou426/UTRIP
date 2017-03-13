c======================================================================|
      subroutine getfbx(fx,fy,ez,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ez(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      do j=1,jx
      do i=1,ix
         fx(i,j)= 0.d0
         fy(i,j)= ez(i,j)
      enddo
      enddo

      return
      end
