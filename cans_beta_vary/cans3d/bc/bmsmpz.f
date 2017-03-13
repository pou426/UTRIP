c======================================================================|
      subroutine bmsmpz(mbnd,margin,cmat,ix,jx,kx)
c======================================================================|
c
c NAME  bmsmpz
c
c PURPOSE
c    apply symmetric boundary condition for the coefficient matrix 
c    of the descretized heat conduction equation.
c    The symmetric point is on the grid point.
c    The values out of the boundary have same sign.
c
c INPUTS & OUTPUTS
c    cmat(ix,jx,kx,7): [double] coefficient matrix of heat conduction eq
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
      dimension cmat(ix,jx,kx,7)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then
        kbnd=1+margin
        do j=1,jx
        do i=1,ix
          cmat(i,j,kbnd,7) = cmat(i,j,kbnd,6)+cmat(i,j,kbnd,7)
          cmat(i,j,kbnd,6) = 0.d0
        enddo
        enddo
      else
        kbnd=kx-margin
        do j=1,jx
        do i=1,ix
          cmat(i,j,kbnd,6) = cmat(i,j,kbnd,6)+cmat(i,j,kbnd,7)
          cmat(i,j,kbnd,7) = 0.d0
        enddo
        enddo
      endif

      return
      end
