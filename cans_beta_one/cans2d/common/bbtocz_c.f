c======================================================================|
      subroutine bbtocz_c(cz,by,dx,x,ix,jx)
c======================================================================|
c
c NAME  bbtocz_c
c
c PURPOSE
c    calculate current density
c        * z-component in cylindrical coordinate
c
c OUTPUTS
c    cz(ix,jx): [double] current density
c
c INPUTS
c    by(ix,jx): [double] magnetic field
c    dx(ix) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    bug-fixed 2006-2-28 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension x(ix)
      dimension dx(ix)
      dimension by(ix,jx)
      dimension cz(ix,jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cz(i,j) =  (x(i+1)*by(i+1,j)-x(i-1)*by(i-1,j))/dx(i)/2./x(i)
      enddo
      enddo

      return
      end
