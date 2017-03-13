c======================================================================|
      subroutine uwd_b(u,dt,dx,ix)
c======================================================================|
c
c NAME  uwd_b
c
c PURPOSE
c    solve Burgers eq. by upwind method
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

      dimension dx(ix)
      dimension u(ix)
      dimension fx(ix)
      dimension cx(ix),f0x(ix)
c----------------------------------------------------------------------|

      do i=1,ix-1
         f0x(i)= u(i)*u(i)/2
      enddo

      do i=1,ix-1
         cx(i)= (u(i+1)+u(i))/2
      enddo

      do i=1,ix-1
         fx(i)= ( ( f0x(i+1)+f0x(i) ) - abs(cx(i))*( u(i+1)-u(i) ) )/2
      enddo


      do i=2,ix-1
         u(i)= u(i) -dt/dx(i)*(fx(i)-fx(i-1))
      enddo

      return
      end
