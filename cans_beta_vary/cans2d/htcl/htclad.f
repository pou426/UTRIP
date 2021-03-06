c======================================================================|
      subroutine htclad(te,dt,gm,cool,htst,ix,jx)
c======================================================================|
c    
c NAME  htclad
c    
c PURPOSE
c    add radiative cooling & static heating term
c    
c INPUTS & OUTPUTS
c    te(ix,jx): [double] temperature
c    
c OUTPUTS
c    None
c
c INPUTS
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    cool(ix,jx): [double] radiative cooling term
c    htst(ix,jx): [double] static heating
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension te(ix,jx)
      dimension cool(ix,jx)
      dimension htst(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         dte=htst(i,j)-cool(i,j)
         te(i,j)=te(i,j)+dt*(gm-1)*dte
      enddo
      enddo

      return
      end
