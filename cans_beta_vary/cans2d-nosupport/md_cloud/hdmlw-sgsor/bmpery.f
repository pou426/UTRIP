c======================================================================|
      subroutine bmpery(margin0,margin1,cmat,ix,jx)
c======================================================================|
c
c NAME  bmsmpy
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

        jbnd0=1+margin0
        jbnd1=jx-margin1

        do j=1,margin0
        do i=1,ix
          cmat(i,jbnd0-j,1) = cmat(i,jbnd1+1-j,1)
          cmat(i,jbnd0-j,2) = cmat(i,jbnd1+1-j,2)
          cmat(i,jbnd0-j,3) = cmat(i,jbnd1+1-j,3)
          cmat(i,jbnd0-j,4) = cmat(i,jbnd1+1-j,4)
          cmat(i,jbnd0-j,5) = cmat(i,jbnd1+1-j,5)
        enddo
        enddo

        do j=1,margin1
        do i=1,ix
          cmat(i,jbnd1+j,1) = cmat(i,jbnd0-1+j,1)
          cmat(i,jbnd1+j,2) = cmat(i,jbnd0-1+j,2)
          cmat(i,jbnd1+j,3) = cmat(i,jbnd0-1+j,3)
          cmat(i,jbnd1+j,4) = cmat(i,jbnd0-1+j,4)
          cmat(i,jbnd1+j,5) = cmat(i,jbnd0-1+j,5)
        enddo
        enddo

      return
      end
