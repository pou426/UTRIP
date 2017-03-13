c======================================================================|
      subroutine bmsmpx(mbnd,margin,cmat,ix,jx,kx)
c======================================================================|
c
c NAME  bmsmpx
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
        ibnd=1+margin
        do k=1,kx
        do j=1,jx
          cmat(ibnd,j,k,3) = cmat(ibnd,j,k,2)+cmat(ibnd,j,k,3)
          cmat(ibnd,j,k,2) = 0.d0
        enddo
        enddo
      else
        ibnd=ix-margin
        do k=1,kx
        do j=1,jx
          cmat(ibnd,j,k,2) = cmat(ibnd,j,k,2)+cmat(ibnd,j,k,3)
          cmat(ibnd,j,k,3) = 0.d0
        enddo
        enddo
      endif

      return
      end
