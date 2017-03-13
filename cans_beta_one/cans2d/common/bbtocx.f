c======================================================================|
      subroutine bbtocx(cx,bz,dy,ix,jx)
c======================================================================|
c
c NAME  bbtocx
c
c PURPOSE
c    calculate current density
c        * x-component
c
c OUTPUTS
c    cx(ix,jx): [double] current density
c
c INPUTS
c    bz(ix,jx): [double] magnetic field
c    dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dy(jx)
      dimension bz(ix,jx)
      dimension cx(ix,jx)
c----------------------------------------------------------------------|

         do j=2,jx-1
         do i=2,ix-1
           cx(i,j) =  (bz(i,j+1)-bz(i,j-1))/dy(j)/2.
         enddo
         enddo

          return
          end
