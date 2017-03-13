c======================================================================|
      subroutine prtote(te,ro,pr,gm,ix,jx,kx)
c======================================================================|
c 
c NAME  prtote
c    
c PURPOSE
c    calculation of temperature
c    
c OUTPUTS
c    te(ix,jx,kx): [double] temperature
c    
c INPUTS
c    ro(ix,jx,kx): [double] density
c    pr(ix,jx,kx): [double] pressure
c    gm: [double] polytropic index gamma
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c 
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx,kx)
      dimension pr(ix,jx,kx)
      dimension te(ix,jx,kx)
c----------------------------------------------------------------------|

         do k=1,kx
         do j=1,jx
         do i=1,ix
            te(i,j,k)=gm*pr(i,j,k)/ro(i,j,k)
         enddo
         enddo
         enddo

         return
         end
