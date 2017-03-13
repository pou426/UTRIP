c======================================================================|
      subroutine htcl(ro,pr,dt,gm,cool0,tecl0,decl0,htst,ix)
c======================================================================|
c
c NAME  htcl
c
c PURPOSE
c    calculate radiative cooling & static heating
c
c INPUTS & OUTPUTS
c    pr(ix): [double] pressure
c
c OUTPUTS
c    None
c
c INPUTS
c    ro(ix): [double] density
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c                cooling time is tau_cool=tecl0/cool0
c    decl0: [double] threshold for optically-thick
c    htst(ix): [double] static heating
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
      dimension htst(ix)
      dimension cool(ix)
c----------------------------------------------------------------------|

        call prtote(te,ro,pr,gm,ix)

        call cooldef(cool,cool0,tecl0,decl0,ro,te,ix)
        call htclad(te,dt,gm,cool,htst,ix)

        call tetopr(ro,pr,te,gm,ix)

      return
      end
