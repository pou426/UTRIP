c======================================================================|
      subroutine bbtocy(cy,bz,dx,ix,jx)
c======================================================================|
c
c NAME  bbtocy
c
c PURPOSE
c    calculate current density
c        * y-component
c
c OUTPUTS
c    cy(ix,jx): [double] current density
c
c INPUTS
c    bz(ix,jx): [double] magnetic field
c    dx(ix) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix)
      dimension bz(ix,jx)
      dimension cy(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cy(i,j) = -(bz(i+1,j)-bz(i-1,j))/dx(i)/2.
      enddo
      enddo

      return
      end
