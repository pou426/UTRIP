c======================================================================|
      subroutine tetopr(ro,pr,te,gm,ix)
c======================================================================|
c
c NAME  tetopr
c
c PURPOSE
c    calculation of pressure from temperature
c
c OUTPUTS
c    pr(ix): [double] pressure
c
c INPUTS
c    te(ix): [double] temperature
c    ro(ix): [double] density
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
        pr(i)=ro(i)*te(i)/gm
      enddo

      return
      end
