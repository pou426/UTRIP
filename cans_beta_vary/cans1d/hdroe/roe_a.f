c======================================================================|
      subroutine roe_a(ro,dt,vx,dx,ix)
c======================================================================|
c
c NAME  roe_a
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * simple advection
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vx(ix), vxm(ix) : [double] velocity along the x-cordinate
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)

      dimension ro(ix),vx(ix)

      dimension fro(ix)
      dimension roh(ix)
      dimension row(ix,2),vxw(ix,2)

c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do i=1,ix-1
         row(i,1)=ro(i)
         vxw(i,1)=vx(i)
         row(i,2)=ro(i+1)
         vxw(i,2)=vx(i+1)
      enddo

      call roeflux_a(fro,row,vxw,ix)

      do i=2,ix-1
         roh(i)=ro(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
      enddo

c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)

      call roeflux_a(fro,row,vxw,ix)

      do i=3,ix-2
         ro(i)=ro(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
      enddo
c----------------------------------------------------------------------|
      return
      end
