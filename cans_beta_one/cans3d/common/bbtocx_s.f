c======================================================================|
      subroutine bbtocx_s(cx,by,bz,dy,dz,x,y,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocx_s
c
c PURPOSE
c    calculate current density
c        * r-component in spherical coordinate
c
c OUTPUTS
c    cx(ix,jx,kx): [double] current density
c
c INPUTS
c    by(ix,jx,kx): [double] magnetic field
c    dz(kx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dy(jx),dz(kx)
      dimension by(ix,jx,kx),bz(ix,jx,kx)
      dimension cx(ix,jx,kx)
      dimension x(ix),y(jx)
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cx(i,j,k) = (
     &      (bz(i,j+1,k)*sin(y(j+1))-bz(i,j-1,k)*sin(y(j-1)))/dy(j)/2.
     &      -(by(i,j,k+1)-by(i,j,k-1))/dz(k)/2.
     &                  ) /x(i)/sin(y(j))
         enddo
         enddo
         enddo

          return
          end
