c======================================================================|
      subroutine cipdxsrc(dadx,da,dah,u,dt,dx,ix)
c======================================================================|
c
c NAME  cipdxsrc
c
c PURPOSE
c    advance no-advective (source-term) phase of CIP method for
c    physical variable gradients
c
c INPUTS & OUTPUTS
c    dadx(ix): [double] physical variable gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    da(ix): [double] physical variable
c    dah(ix): [double] physical variable after non-advective phase
c    u(ix): [double] advection velocity
c    dt: [double] delta time
c    dx(ix) : [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dadx(ix),u(ix),dx(ix)
      dimension da(ix),dah(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1
        dadx(i)=dadx(i)
     &           +((dah(i+1)-da(i+1))-(dah(i-1)-da(i-1)))/2/dx(i)
     &           -dadx(i)*(u(i+1)-u(i-1))*dt/dx(i)/2 
      enddo

      return
      end
