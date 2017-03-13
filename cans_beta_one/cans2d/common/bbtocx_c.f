c======================================================================|
      subroutine bbtocx_c(cx,by,dz,ix,jx)
c======================================================================|
c
c NAME  bbtocx_c
c
c PURPOSE
c    calculate current density
c        * r-component in cylindrical coordinate
c
c OUTPUTS
c    cx(ix,jx): [double] current density
c
c INPUTS
c    by(ix,jx): [double] magnetic field
c    dz(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dz(jx)
      dimension by(ix,jx)
      dimension cx(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cx(i,j) = -(by(i,j+1)-by(i,j-1))/dz(j)/2.
      enddo
      enddo

      return
      end
