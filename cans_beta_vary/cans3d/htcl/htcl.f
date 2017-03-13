c======================================================================|
      subroutine htcl(ro,pr,dt,gm,cool0,tecl0,decl0,htst,ix,jx,kx)
c======================================================================|
c
c NAME  htcl
c
c PURPOSE
c    calculate radiative cooling & static heating 
c
c INPUTS & OUTPUTS
c    pr(ix,jx): [double] pressure
c
c OUTPUTS
c    None
c    
c INPUTS         
c    ro(ix,jx): [double] density
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c                cooling time is tau_cool=tecl0/cool0
c    decl0: [double] threshold for optically-thick
c    htst(ix,jx): [double] static heating
c    ix,jx: [integer] dimension size
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx)
      dimension pr(ix,jx,kx)
      dimension te(ix,jx,kx)
      dimension htst(ix,jx,kx)
      dimension cool(ix,jx,kx)
c----------------------------------------------------------------------|

        call prtote(te,ro,pr,gm,ix,jx,kx)

        call cooldef(cool,cool0,tecl0,decl0,ro,te,ix,jx,kx)
        call htclad(te,dt,gm,cool,htst,ix,jx,kx)

        call tetopr(ro,pr,te,gm,ix,jx,kx)

      return
      end
