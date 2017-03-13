c======================================================================|
      subroutine htclprm(cool0,tecl0,decl0,gm,tenml,ronml,rlnml)
c======================================================================|
c
c NAME  htclprm
c
c PURPOSE
c    calculate non-dimensional conduction coefficient from
c    dimensional parameters
c 
c OUTPUTS
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c               cooling time is tau_cool=tecl0/cool0 @ (te=tecl0, de=1)
c    decl0: [double] threshold for optically-thick
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

      rkb=1.38e-16
      rmm=1.67e-24

      venml=sqrt(gm*rkb/rmm*tenml)
      tmnml=rlnml/venml
      rknml=rkb*tenml/(ronml*tmnml)

c  cooling function
c   cooling time is : tau_cool=tecl0/cool0 @ (te=tecl0, de=1)

      cool0=8.0e-22/rknml
      tecl0=2.d5/tenml

c  critical density for optically-thick condition

      decl0=1.d12/ronml


      return
      end
