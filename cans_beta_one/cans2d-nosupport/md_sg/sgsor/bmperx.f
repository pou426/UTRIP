c======================================================================|
      subroutine bmperx(margin0,margin1,cmat,ix,jx)
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
c    cmat(ix,jx,5): [double] coefficient matrix of heat conduction eq
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
      dimension cmat(ix,jx,5)
c----------------------------------------------------------------------|

        ibnd0=1+margin0
        ibnd1=ix-margin1

        do i=1,margin0
        do j=1,jx
          cmat(ibnd0-i,j,1) = cmat(ibnd1+1-i,j,1)
          cmat(ibnd0-i,j,2) = cmat(ibnd1+1-i,j,2)
          cmat(ibnd0-i,j,3) = cmat(ibnd1+1-i,j,3)
          cmat(ibnd0-i,j,4) = cmat(ibnd1+1-i,j,4)
          cmat(ibnd0-i,j,5) = cmat(ibnd1+1-i,j,5)
        enddo
        enddo

        do i=1,margin1
        do j=1,jx
          cmat(ibnd1+i,j,1) = cmat(ibnd0-1+i,j,1)
          cmat(ibnd1+i,j,2) = cmat(ibnd0-1+i,j,2)
          cmat(ibnd1+i,j,3) = cmat(ibnd0-1+i,j,3)
          cmat(ibnd1+i,j,4) = cmat(ibnd0-1+i,j,4)
          cmat(ibnd1+i,j,5) = cmat(ibnd0-1+i,j,5)
        enddo
        enddo

      return
      end
