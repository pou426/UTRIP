c======================================================================|
      subroutine getsaz(ss,ez,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ez(ix,jx)
      dimension ss(ix,jx)

      do j=1,jx
      do i=1,ix
         ss(i,j)= -ez(i,j)
      enddo
      enddo

      return
      end
