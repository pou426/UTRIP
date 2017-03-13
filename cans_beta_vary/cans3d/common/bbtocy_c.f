c======================================================================|
      subroutine bbtocy_c(cy,bz,bx,dz,dx,ix,jx,kx)
c======================================================================|
c
c NAME  bbtocy_c
c
c PURPOSE
c    calculate current density
c        * phi-component in cylindrical coordinate
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
      dimension bz(ix,jx,kx),bx(ix,jx,kx)
      dimension cy(ix,jx,kx)
c----------------------------------------------------------------------|

         do k=2,kx-1
         do j=2,jx-1
         do i=2,ix-1
           cy(i,j,k) =  (bx(i,j,k+1)-bx(i,j,k-1))/dz(k)/2.
     &                 -(bz(i+1,j,k)-bz(i-1,j,k))/dx(i)/2.
         enddo
         enddo
         enddo

          return
          end
