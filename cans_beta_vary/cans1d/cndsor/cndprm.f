c======================================================================|
      subroutine cndprm(rkap0,gm,tenml,ronml,rlnml)
c======================================================================|
c
c NAME  cndprm
c
c PURPOSE
c    calculate non-dimensional conduction coefficient from
c    dimensional parameters  
c
c OUTPUTS
c    rkap0: [double] constant part of heat conduction coefficient
c                conduction time is tau_cnd=1/rkap0 @ (te=1,de=1)
c
c INPUTS
c    gm: [double] polytropic index gamma
c    tenml: [double] normalization unit for temperature in Kelvin
c    ronml: [double] normalization unit for density in cm^-3
c    rlnml: [double] normalization unit for length in cm
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      rkb=1.38d-16
      rmm=1.67d-24
      rka=1.d-6

      venml=sqrt(gm*rkb/rmm*tenml)
      rkpnml=rmm*ronml*venml*rlnml*rkb/rmm
      rkap0=rka/rkpnml*sqrt(tenml)**5


      return
      end
