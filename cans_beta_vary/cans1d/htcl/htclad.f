c======================================================================|
      subroutine htclad(te,dt,gm,cool,htst,ix)
c======================================================================|
c
c NAME  htclad
c
c PURPOSE
c    add radiative cooling & static heating term
c
c INPUTS & OUTPUTS
c    te(ix): [double] temperature
c
c OUTPUTS
c    None
c
c INPUTS
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    cool(ix): [double] radiative cooling term
c    htst(ix): [double] static heating
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension te(ix)
      dimension cool(ix)
      dimension htst(ix)
c----------------------------------------------------------------------|

      do i=1,ix
         dte=htst(i)-cool(i)
         te(i)=te(i)+dt*(gm-1)*dte
      enddo


      return
      end
