c======================================================================|
      subroutine daupdate(u,un,ix,jx,kx)
c======================================================================|

      implicit double precision (a-h,o-z)
      dimension u(ix,jx,kx), un(ix,jx,kx)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         u(i,j,k)= un(i,j,k)
      enddo
      enddo
      enddo

      return
      end
