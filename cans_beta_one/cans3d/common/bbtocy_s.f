c======================================================================|
      subroutine bbtocy_s(cy,bz,bx,dz,dx,x,y,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocy_s
c
c PURPOSE
c    calculate current density
c        * theta-component in spherical coordinate
c
c OUTPUTS
c    cy(ix,jx,kx): [double] current density
c
c INPUTS
c    bx(ix,jx,kx): [double] magnetic field
c    bz(ix,jx,kx): [double] magnetic field
c    dx(ix),dz(jx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dz(kx),dx(ix)
      dimension x(ix),y(jx)
      dimension bz(ix,jx,kx),bx(ix,jx,kx)
      dimension cy(ix,jx,kx)
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cy(i,j,k) = (bx(i,j,k+1)-bx(i,j,k-1))/dz(k)/2./x(i)/sin(y(j))
     &           -(bz(i+1,j,k)*x(i+1)-bz(i-1,j,k)*x(i-1))/dx(i)/2./x(i)
         enddo
         enddo
         enddo

          return
          end
