c======================================================================|
      subroutine bbtocy_s(cy,bz,dx,x,y,ix,jx)
c======================================================================|
c
c NAME  bbtocy_s
c
c PURPOSE
c    calculate current density
c        * phi-component in spherical coordinate
c
c OUTPUTS
c    cy(ix,jx): [double] current density
c
c INPUTS
c    bz(ix,jx): [double] magnetic field
c    x(ix),y(jx) : [double] coordinate
c    dx(ix) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx)
      dimension dx(ix)
      dimension bz(ix,jx)
      dimension cy(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cy(i,j) = -(x(i+1)*sin(y(j))*bz(i+1,j)
     &             -x(i-1)*sin(y(j))*bz(i-1,j))/dx(i)/2.
     &          /x(i)/sin(y(j))
      enddo
      enddo

      return
      end
