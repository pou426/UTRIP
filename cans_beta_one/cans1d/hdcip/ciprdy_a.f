c======================================================================|
      subroutine ciprdy_a(ro,rodx,dx,ix)
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
c    rodx(ix): [double] density gradient
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
c    dx(ix) : [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)
      dimension ro(ix),rodx(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1
        rodx(i)=(ro(i+1)-ro(i-1))/dx(i)/2
      enddo


      return
      end
