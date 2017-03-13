c======================================================================|
      subroutine bdcndx(mbnd,margin,qq,dqq,ix,jx,kx)
c======================================================================|
c
c NAME  bdcndx
c
c PURPOSE
c    apply constant-gradient boundary condition
c    values are extended to have constant gradient
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
      dimension qq(ix,jx,kx),dqq(margin,jx,kx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        do k=1,kx 
        do j=1,jx 
        do i=1,margin 
          qq(ibnd-i,j,k) = qq(ibnd,j,k)-dqq(ibnd-i,j,k)
        enddo
        enddo
        enddo
      else
        ibnd=ix-margin
        do k=1,kx 
        do j=1,jx 
        do i=1,margin
          qq(ibnd+i,j,k) = qq(ibnd,j,k)+dqq(ibnd+i,j,k)
        enddo
        enddo
        enddo
      endif
      
      return
      end
