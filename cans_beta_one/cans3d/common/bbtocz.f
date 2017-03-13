c======================================================================|
      subroutine bbtocz(cz,bx,by,dx,dy,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocz
c
c PURPOSE
c    calculate current density
c        * z-component
c
c OUTPUTS
c    cz(ix,jx,kx): [double] current density
c
c INPUTS
c    bx(ix,jx,kx): [double] magnetic field
c    by(ix,jx,kx): [double] magnetic field
c    dx(ix),dy(jx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension bx(ix,jx,kx),by(ix,jx,kx)
      dimension cz(ix,jx,kx)
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cz(i,j,k) =  (by(i+1,j,k)-by(i-1,j,k))/dx(i)/2.
     &                 -(bx(i,j+1,k)-bx(i,j-1,k))/dy(j)/2.
         enddo
         enddo
         enddo

          return
          end
