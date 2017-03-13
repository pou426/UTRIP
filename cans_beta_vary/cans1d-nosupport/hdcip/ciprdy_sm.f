c======================================================================|
      subroutine ciprdy_sm(de,ei,rxm,ry,dedx,eidx,rxdxm,rydx,dx,dxm,ix)
c======================================================================|
c
c NAME  ciprdy_h
c
c PURPOSE
c    derive temperature: derive velocity between grid points
c    derive gradient of physical variables
c        * hydrodynamics
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    te(ix): [double] temperature
c    vxm(ix): [double] velocity
c    rodx(ix): [double] density gradient
c    tedx(ix): [double] temperature gradient
c    vxdxm(ix): [double] velocity gradient
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity
c    gm: [double] polytropic index gamma
c    dx(ix),dxm(ix) : [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)

      dimension de(ix),ei(ix),rxm(ix),ry(ix)
      dimension dedx(ix),eidx(ix),rxdxm(ix),rydx(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1
        dedx(i)=(de(i+1)-de(i-1))/dx(i)/2
        eidx(i)=(ei(i+1)-ei(i-1))/dx(i)/2
        rydx(i)=(ry(i+1)-ry(i-1))/dx(i)/2
        rxdxm(i)=(rxm(i+1)-rxm(i-1))/dxm(i)/2
      enddo


      return
      end
