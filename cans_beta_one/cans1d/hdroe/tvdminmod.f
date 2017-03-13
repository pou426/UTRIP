c======================================================================|
      subroutine tvdminmod(da,daw,ix)
c======================================================================|
c
c NAME  tvdminmod
c
c PURPOSE
c    Interporate the physical variables based on MUSCL
c    using 'min-mod' function as a limitter
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    daw(ix,2): [double] variable at cell boundary
c
c INPUTS
c    da(ix): [double] physical variable
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension da(ix)
      dimension daw(ix,2)

c----------------------------------------------------------------------|
c     define limiter functions

      flmt(a,b)=max(0.0d0,min(b*sign(1.0d0,a),abs(a)))*sign(1.0d0,a)

c----------------------------------------------------------------------|

      do i=2,ix-2
         daw(i,1)=da(i)+0.5*flmt(da(i+1)-da(i),da(i)-da(i-1))
         daw(i,2)=da(i+1)-0.5*flmt(da(i+1)-da(i),da(i+2)-da(i+1))
      enddo

      return
      end
