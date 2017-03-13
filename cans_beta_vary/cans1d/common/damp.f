c======================================================================|
      subroutine damp(da,xdamp,wdamp,damp0,x,ix)
c======================================================================|
c
c NAME  damp
c
c PURPOSE
c    set damping zone to reduce the wave amplitude.
c
c INPUTS & OUTPUTS
c    da(ix): [double] variable
c
c OUTPUTS
c    None
c
c INPUTS
c    xdamp: [double] central coordinate of damping zone
c    wdamp: [double] width of damping zone
c    damp0: [double] damping strength
c    x(ix): [double] coordinate
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension da(ix)
      dimension x(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1
        cdamp=damp0*exp(-(x(i)-xdamp)**2/wdamp**2)
        da(i)=da(i)+cdamp*(da(i-1)-2*da(i)+da(i+1))
      enddo


      return
      end
