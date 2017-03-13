c======================================================================|
      subroutine tetopr(ro,pr,te,gm,ix,jx,kx)
c======================================================================|
c 
c NAME  tetopr
c    
c PURPOSE
c    calculation of pressure from temperature
c    
c OUTPUTS
c    pr(ix,jx,kx): [double] pressure
c    
c INPUTS
c    te(ix,jx,kx): [double] temperature
c    ro(ix,jx,kx): [double] density
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
        pr(i,j,k)=ro(i,j,k)*te(i,j,k)/gm
      enddo
      enddo
      enddo

      return
      end
