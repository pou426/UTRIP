c======================================================================|
      subroutine bbtocz_s(cz,bx,by,dx,dy,x,ix,jx)
c======================================================================|
c
c NAME  bbtocz_s
c
c PURPOSE
c    calculate current density
c        * theta-component in spherical coordinate
c
c OUTPUTS
c    cz(ix,jx): [double] current density
c
c INPUTS
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    x(ix) : [double] coordinate
c    dx(ix),dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension x(ix)
      dimension dx(ix),dy(jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension cz(ix,jx)

c----------------------------------------------------------------------|

         do j=2,jx-1
         do i=2,ix-1
           cz(i,j) =( (x(i+1)*by(i+1,j)-x(i-1)*by(i-1,j))/dx(i)/2.
     &               -(bx(i,j+1)-bx(i,j-1))/dy(j)/2)
     &               /x(i)
         enddo
         enddo

          return
          end
