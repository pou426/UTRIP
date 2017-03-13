c======================================================================|
      subroutine daupdate(u,un,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)
      dimension u(ix,jx), un(ix,jx)

      do j=1,jx
      do i=1,ix
         u(i,j)= un(i,j)
      enddo
      enddo

      return
      end
