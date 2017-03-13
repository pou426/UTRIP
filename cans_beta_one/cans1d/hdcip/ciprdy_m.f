c======================================================================|
      subroutine ciprdy_m(te,vxm,rodx,tedx,vxdxm,vydx,ro,pr,vx,vy
     &    ,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  ciprdy_m
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
c    vydx(ix): [double] velocity gradient
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity
c    vy(ix): [double] velocity
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

      dimension ro(ix),pr(ix),vx(ix),vy(ix)
      dimension te(ix),vxm(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix),vydx(ix)
c----------------------------------------------------------------------|

      do i=1,ix
         te(i)=gm*pr(i)/ro(i)
      enddo

      do i=1,ix-1
        vxm(i)=(vx(i)+vx(i+1))/2
      enddo

      do i=2,ix-1
        rodx(i)=(ro(i+1)-ro(i-1))/dx(i)/2
        tedx(i)=(te(i+1)-te(i-1))/dx(i)/2
        vydx(i)=(vy(i+1)-vy(i-1))/dx(i)/2
        vxdxm(i)=(vxm(i+1)-vxm(i-1))/dxm(i)/2
      enddo


      return
      end
