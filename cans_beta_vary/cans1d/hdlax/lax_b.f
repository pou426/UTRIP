c======================================================================|
      subroutine lax_b(u,dt,dx,dxm,ix)
c======================================================================|
c
c NAME  lax_b
c
c PURPOSE
c    solve Burgers eq. by Lax-Friedrich method
c
c INPUTS & OUTPUTS
c    u(ix): [double] 
c
c OUTPUTS
c    None
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
      dimension fx(ix)
c----------------------------------------------------------------------|

      do i=1,ix
         fx(i)= u(i)*u(i)/2
      enddo
      call lax(un,u,dt,fx,dx,dxm,ix)
      do i=2,ix-1
         u(i) = un(i)
      enddo

      return
      end
