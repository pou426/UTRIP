c======================================================================|
      subroutine bbtocx(cx,by,bz,dy,dz,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocx
c
c PURPOSE
c    calculate current density
c        * x-component
c
c OUTPUTS
c    cx(ix,jx,kx): [double] current density
c
c INPUTS
c    bz(ix,jx,kx): [double] magnetic field
c    dy(jx) : [double] grid spacing
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
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cx(i,j,k) =  (bz(i,j+1,k)-bz(i,j-1,k))/dy(j)/2.
     &                 -(by(i,j,k+1)-by(i,j,k-1))/dz(k)/2.
         enddo
         enddo
         enddo

          return
          end
