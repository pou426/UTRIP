c======================================================================|
      subroutine htstdef(htst,gm,cool0,tecl0,decl0,ro,pr,ix)
c======================================================================|
c
c NAME  htstdef
c
c PURPOSE
c    define static heating term. equal to initial cooling term
c
c OUTPUTS
c    htst(ix): [double] static heating
c
c INPUTS
c    gm: [double] polytropic index gamma
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c               cooling time is tau_cool=tecl0/cool0 @ (te=tecl0, de=1)
c    decl0: [double] threshold for optically-thick
c    ro(ix): [double] density
c    pr(ix): [double] pressure
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
      dimension cool(ix)
      dimension htst(ix)
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix)
      call cooldef(cool,cool0,tecl0,decl0,ro,te,ix)

      do i=1,ix
           htst(i)=cool(i)
      enddo


      return
      end
