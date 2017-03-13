c======================================================================|
      subroutine cip_a(ro,rodx,dt,vx,vxm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_a
c
c PURPOSE
c    solve eqs. by CIP  method with effects of
c        * simple advection
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    rodx(ix): [double] density gradient
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

      dimension dx(ix),dxm(ix)
      dimension vx(ix),vxm(ix)
      dimension ro(ix),rodx(ix)
      dimension roh(ix),rodxh(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix
        roh(i)=ro(i)
        rodxh(i)=rodx(i)
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
      do i=2,ix
        du=vxm(i)-vxm(i-1)
        roh(i)=roh(i)+dt*(-du/dx(i)*ro(i))
      enddo

      call cipdxsrc(rodxh,ro,roh,vx,dt,dx,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      call cipadv(roh,rodxh,vx,0,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do i=1,ix
        ro(i)=roh(i)
        rodx(i)=rodxh(i)
      enddo


      return
      end
