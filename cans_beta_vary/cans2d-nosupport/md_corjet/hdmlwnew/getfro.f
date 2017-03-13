c======================================================================|
      subroutine getfro(fx,fy,ro,vx,vy,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)
      enddo
      enddo

      return
      end
