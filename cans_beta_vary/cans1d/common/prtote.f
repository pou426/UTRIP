c======================================================================|
      subroutine prtote(te,ro,pr,gm,ix)
c======================================================================|
c
c NAME  prtote
c
c PURPOSE
c    calculation of temperature
c
c OUTPUTS
c    te(ix): [double] temperature
c
c INPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension pr(ix)
      dimension te(ix)
c----------------------------------------------------------------------|
         do i=1,ix
            te(i)=gm*pr(i)/ro(i)
         enddo

         return
         end
