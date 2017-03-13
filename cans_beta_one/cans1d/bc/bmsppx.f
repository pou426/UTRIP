c======================================================================|
      subroutine bmsppx(mbnd,margin,cmat,ix)
c======================================================================|
c
c NAME  bdsmnx
c
c PURPOSE
c    apply symmetric boundary condition for the coefficient matrix 
c    of the descretized heat conduction equation.
c    The symmetric point is between the grid point.
c    The values out of the boundary have same sign.
c
c INPUTS & OUTPUTS
c    cmat(ix,3): [double] coefficient matrix of heat conduction eq
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
      dimension cmat(ix,3)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then
        ibnd=1+margin
        cmat(ibnd,1) = cmat(ibnd,1)+cmat(ibnd,2)
        cmat(ibnd,2) = 0.d0
      else
        ibnd=ix-margin
        cmat(ibnd,1) = cmat(ibnd,1)+cmat(ibnd,3)
        cmat(ibnd,3) = 0.d0
      endif

      return
      end
