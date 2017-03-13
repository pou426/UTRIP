c======================================================================|
      subroutine bdspnx(mbnd,margin,qq,ix,jx)
c======================================================================|
c
c NAME  bdspnx
c
c PURPOSE
c    apply symmetric boundary condition. 
c    The symmetric point is between the grid point.
c    The values out of the boundary have opposite sign.
c
c INPUTS & OUTPUTS
c    qq(ix,jx): [double] variable
c
c OUTPUTS
c    None
c
c INPUTS
c    ix,jx: [integer] dimension size
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    mbnd: [integer] If mbnd=0, smaller 'i' side. 
c                    If mbnd=1, larger  'i' side.
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension qq(ix,jx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then
        ibnd=1+margin
        do i=1,margin
        do j=1,jx
          qq(ibnd-i,j) = -qq(ibnd-1+i,j)
        enddo
        enddo
      else
        ibnd=ix-margin
        do i=1,margin
        do j=1,jx
          qq(ibnd+i,j) = -qq(ibnd+1-i,j)
        enddo
        enddo
      endif
      
      return
      end
