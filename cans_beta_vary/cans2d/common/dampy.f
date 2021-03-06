c======================================================================|
      subroutine dampy(da,ydamp,wdamp,damp0,y,ix,jx)
c======================================================================|
c
c NAME  dampy
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
      dimension da(ix,jx)
      dimension y(jx)
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
        cdamp=damp0*exp(-(y(j)-ydamp)**2/wdamp**2)
        da(i,j)=da(i,j)+cdamp*(da(i-1,j)-2*da(i,j)+da(i+1,j))
      enddo
      enddo


      return
      end
