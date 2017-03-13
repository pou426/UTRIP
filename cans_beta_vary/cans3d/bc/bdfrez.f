c======================================================================|
      subroutine bdfrez(mbnd,margin,qq,ix,jx,kx)
c======================================================================|
c
c NAME  bdfrez
c
c PURPOSE
c    apply free boundary condition 
c    
c INPUTS & OUTPUTS
c    qq(ix,jx,kx): [double] variable
c    
c OUTPUTS
c    None
c    
c INPUTS
c    ix,jx,kx: [integer] dimension size
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    mbnd: [integer] If mbnd=0, smaller 'i' side. 
c                    If mbnd=1, larger  'i' side. 
c                    
c HISTORY            
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension qq(ix,jx,kx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        kbnd=1+margin
        do k=1,margin
        do j=1,jx
        do i=1,ix
          qq(i,j,kbnd-k) = qq(i,j,kbnd)
        enddo
        enddo
        enddo
      else
        kbnd=kx-margin
        do k=1,margin
        do j=1,jx
        do i=1,ix
          qq(i,j,kbnd+k) = qq(i,j,kbnd)
        enddo
        enddo
        enddo
      endif
      
      return
      end
