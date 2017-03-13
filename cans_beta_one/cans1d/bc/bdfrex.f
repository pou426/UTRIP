c======================================================================|
      subroutine bdfrex(mbnd,margin,qq,ix)
c======================================================================|
c
c NAME  bdfrex
c
c PURPOSE
c    apply free boundary condition
c
c INPUTS & OUTPUTS
c    qq(ix): [double] variable
c
c OUTPUTS
c    None
c
c INPUTS
c    ix: [integer] dimension size
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    mbnd: [integer] If mbnd=0, smaller 'i' side. 
c                    If mbnd=1, larger  'i' side.
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension qq(ix)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        do i=1,margin 
          qq(ibnd-i) = qq(ibnd)
        enddo
      else
        ibnd=ix-margin
        do i=1,margin
          qq(ibnd+i) = qq(ibnd)
        enddo
      endif
      
      return
      end
