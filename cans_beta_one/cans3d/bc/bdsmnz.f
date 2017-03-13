c======================================================================|
      subroutine bdsmnz(mbnd,margin,qq,ix,jx,kx)
c======================================================================|
c
c NAME  bdsmnz
c
c PURPOSE
c    apply symmetric boundary condition.
c    The symmetric point is on the grid point.
c    The values out of the boundary have opposite sign.
c    The value on the symmetric point is fixed to zero.
c    
c OUTPUTS
c    None
c    
c INPUTS & OUTPUTS
c    qq(ix,jx,kx): [double] variable
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
          qq(i,j,kbnd-k) = -qq(i,j,kbnd+k)
        enddo
        enddo
        enddo
        do j=1,jx
        do i=1,ix
          qq(i,j,kbnd) = 0.
        enddo
        enddo
      else
        kbnd=kx-margin
        do k=1,margin
        do j=1,jx
        do i=1,ix
          qq(i,j,kbnd+k) = -qq(i,j,kbnd-k)
        enddo
        enddo
        enddo
        do j=1,jx
        do i=1,ix
          qq(i,j,kbnd) = 0.
        enddo
        enddo
      endif

      return
      end
