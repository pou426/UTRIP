c======================================================================|
      subroutine htstdef(htst,gm,cool0,tecl0,decl0,ro,pr,ix,jx,kx)
c======================================================================|
c
c NAME  htstdef
c
c PURPOSE
c    define static heating term. equal to initial cooling term
c    
c OUTPUTS
c    htst(ix,jx): [double] static heating
c    
c INPUTS
c    gm: [double] polytropic index gamma
c    cool0: [double] strength of radiative cooling
c    tecl0: [double] temperature for peak of cooling function
c               cooling time is tau_cool=tecl0/cool0 @ (te=tecl0, de=1)
c    decl0: [double] threshold for optically-thick
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
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
      dimension cool(ix,jx,kx)
      dimension htst(ix,jx,kx)
c----------------------------------------------------------------------|

      call prtote(te,ro,pr,gm,ix,jx,kx)
      call cooldef(cool,cool0,tecl0,decl0,ro,te,ix,jx,kx)

      do k=1,kx
      do j=1,jx
      do i=1,ix
           htst(i,j,k)=cool(i,j,k)
      enddo
      enddo
      enddo


      return
      end
