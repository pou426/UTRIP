c======================================================================|
      subroutine bdcnsy(mbnd,margin,qq,q0,ix,jx,kx)
c======================================================================|
c
c NAME  bdcnsy
c
c PURPOSE
c    apply constant-value boundary condition
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
c    q0: [double] constant boundary value to be taken
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension qq(ix,jx,kx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        jbnd=1+margin
        do k=1,kx
        do j=1,margin
        do i=1,ix
          qq(i,jbnd-j,k) = q0
        enddo
        enddo
        enddo
      else
        jbnd=jx-margin
        do k=1,kx
        do j=1,margin
        do i=1,ix
          qq(i,jbnd+j,k) = q0
        enddo
        enddo
        enddo
      endif
      
      return
      end
