c======================================================================|
      subroutine ciprdy_a(ro,rodx,rody,dx,ix,dy,jx)
c======================================================================|
c
c NAME  ciprdy_a
c
c PURPOSE
c    derive gradient of physical variables
c        * simple advection
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    rodx(ix,jx): [double] density gradient
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c
c    ro(ix,jx): [double] density
c    dx(ix,jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)
      dimension ro(ix,jx),rodx(ix,jx),rody(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=2,ix-1
        rodx(i,j)=(ro(i+1,j)-ro(i-1,j))/dx(i)/2
      enddo
      enddo

      do j=2,jx-1
      do i=1,ix
        rody(i,j)=(ro(i,j+1)-ro(i,j-1))/dy(j)/2
      enddo
      enddo


      return
      end
