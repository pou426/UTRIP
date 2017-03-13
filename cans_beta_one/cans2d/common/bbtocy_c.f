c======================================================================|
      subroutine bbtocy_c(cy,bz,bx,dz,dx,ix,jx)
c======================================================================|
c
c NAME  bbtocy_c
c
c PURPOSE
c    calculate current density
c        * phi-component in cylindrical coordinate
c
c OUTPUTS
c    cy(ix,jx): [double] current density
c
c INPUTS
c    bx(ix,jx): [double] magnetic field
c    bz(ix,jx): [double] magnetic field
c    dx(ix),dz(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dz(jx)
      dimension bx(ix,jx),bz(ix,jx)
      dimension cy(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cy(i,j) =  (bx(i,j+1)-bx(i,j-1))/dz(j)/2.
     &            -(bz(i+1,j)-bz(i-1,j))/dx(i)/2.
      enddo
      enddo

      return
      end
