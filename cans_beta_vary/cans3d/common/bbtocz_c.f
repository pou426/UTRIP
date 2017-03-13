c======================================================================|
      subroutine bbtocz_c(cz,bx,by,dx,dy,x,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocz_c
c
c PURPOSE
c    calculate current density
c        * z-component in cylindrical coordinate
c
c OUTPUTS
c    cz(ix,jx,kx): [double] current density
c
c INPUTS
c    by(ix,jx,kx): [double] magnetic field
c    dx(ix) : [double] grid spacing
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
      dimension x(ix)
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cz(i,j,k) = 
     &           (x(i+1)*by(i+1,j,k)-x(i-1)*by(i-1,j,k))/dx(i)/2./x(i)
     &                 -(bx(i,j+1,k)-bx(i,j-1,k))/dy(j)/2./x(i)
         enddo
         enddo
         enddo

          return
          end
