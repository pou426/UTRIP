c======================================================================|
      subroutine bdfrdx(mbnd,margin,qq,dxm,ix,jx,kx)
c======================================================================|
c
c NAME  bdfrdx
c
c PURPOSE
c    apply free boundary condition
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
      dimension qq(ix,jx,kx),dxm(ix)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        do k=1,kx 
        do i=1,margin 
        do j=1,jx 
          dqq=(qq(ibnd+1,j,k)-qq(ibnd,j,k))/dxm(ibnd)
          qq(ibnd-i,j,k) = qq(ibnd-i+1,j,k)-dqq*dxm(ibnd-i)
        enddo
        enddo
        enddo
      else
        ibnd=ix-margin
        do k=1,kx 
        do i=1,margin
        do j=1,jx 
          dqq=(qq(ibnd,j,k)-qq(ibnd-1,j,k))/dxm(ibnd-1)
          qq(ibnd+i,j,k) = qq(ibnd+i-1,j,k)+dqq*dxm(ibnd+i-1)
        enddo
        enddo
        enddo
      endif
      
      return
      end
