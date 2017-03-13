c======================================================================|
      subroutine lax(un,u,dt,fx,dx,dxm,ix)
c======================================================================|
c
c NAME  lax
c
c PURPOSE
c    solve Lax-Friedrich method
c
c INPUTS & OUTPUTS
c    u(ix): [double] 
c
c OUTPUTS
c    un(ix): [double] 
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2006-8-2 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension u(ix),un(ix)
      dimension fx(ix),fnx(ix)
c----------------------------------------------------------------------|

      do i=1,ix-1
         fnx(i)= ( ( fx(i+1)+fx(i) ) - dxm(i)/dt*( u(i+1)-u(i) ) )/2
      enddo


      do i=2,ix-1
         un(i)= u(i) -dt/dx(i)*(fnx(i)-fnx(i-1))
      enddo

      return
      end
