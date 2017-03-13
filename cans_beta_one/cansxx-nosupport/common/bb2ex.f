c======================================================================|
      subroutine bb2ex(ex,vy,vz,by,bz,ix,jx,kx)
c======================================================================|
      implicit double precision (a-h,o-z)

      dimension ex(ix,jx,kx)
      dimension by(ix,jx,kx),bz(ix,jx,kx)
      dimension vy(ix,jx,kx),vz(ix,jx,kx)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ex(i,j,k) = -vy(i,j,k)*bz(i,j,k)+vz(i,j,k)*by(i,j,k)
      enddo
      enddo
      enddo


      return
      end
