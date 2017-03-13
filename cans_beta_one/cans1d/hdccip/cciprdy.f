c======================================================================|
      subroutine cciprdy(rdm,ro,dx,ix)
c======================================================================|
c
c NAME  ccipadvrd
c
c PURPOSE
c    solve eqs. by Conservative CIP
c        * simple advection
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    rodx(ix): [double] density gradient
c    rdm(ix) : [double] integrated mass in a cell
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vx(ix), vxm(ix) : [double] velocity along the x-cordinate
c    dx(ix), dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix)
      dimension ro(ix)
      dimension rdm(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix-1
        rdm(i)=(ro(i)*dx(i)+ro(i+1)*dx(i+1))/2.d0
      enddo


      return
      end
