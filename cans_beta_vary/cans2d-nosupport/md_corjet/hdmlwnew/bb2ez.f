c======================================================================|
      subroutine bb2ez(ez,vx,vy,bx,by,ix,jx)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ez(ix,jx)
      dimension vx(ix,jx),vy(ix,jx),bx(ix,jx),by(ix,jx)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         ez(i,j) = -vx(i,j)*by(i,j)+vy(i,j)*bx(i,j)
      enddo
      enddo


      return
      end
