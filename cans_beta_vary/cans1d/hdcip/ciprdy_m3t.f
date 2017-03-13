c======================================================================|
      subroutine ciprdy_m3t(vxm,rodx,vxdxm,vydx,vzdx
     &    ,ro,vx,vy,vz,dx,dxm,ix)
c======================================================================|
c
c NAME  ciprdy_m3
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
c    vxm(ix): [double] velocity
c    rodx(ix): [double] density gradient
c    vxdxm(ix): [double] velocity gradient
c    vydx(ix): [double] velocity gradient
c    vzdx(ix): [double] velocity gradient
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
c    vx(ix): [double] velocity
c    vy(ix): [double] velocity
c    vz(ix): [double] velocity
c    dx(ix),dxm(ix) : [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension ro(ix),vx(ix),vy(ix),vz(ix)
      dimension vxm(ix)
      dimension rodx(ix),vxdxm(ix),vydx(ix),vzdx(ix)
c----------------------------------------------------------------------|

      do i=1,ix-1
        vxm(i)=(vx(i)+vx(i+1))/2
      enddo

      do i=2,ix-1
        rodx(i)=(ro(i+1)-ro(i-1))/dx(i)/2
        vydx(i)=(vy(i+1)-vy(i-1))/dx(i)/2
        vzdx(i)=(vz(i+1)-vz(i-1))/dx(i)/2
        vxdxm(i)=(vxm(i+1)-vxm(i-1))/dxm(i)/2
      enddo


      return
      end
