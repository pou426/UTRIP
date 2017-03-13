c======================================================================|
      subroutine bdfrdx(mbnd,margin,qq,dxm,ix,jx)
c======================================================================|
c
c NAME  bdfrdx
c
c PURPOSE
c    apply free boundary condition
c    values are extended to have constant gradient
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
      dimension qq(ix,jx),dxm(ix)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        do i=1,margin 
        do j=1,jx 
          dqq=(qq(ibnd+1,j)-qq(ibnd,j))/dxm(ibnd)
          qq(ibnd-i,j) = qq(ibnd-i+1,j)-dqq*dxm(ibnd-i)
        enddo
        enddo
      else
        ibnd=ix-margin
        do i=1,margin
        do j=1,jx 
          dqq=(qq(ibnd,j)-qq(ibnd-1,j))/dxm(ibnd-1)
          qq(ibnd+i,j) = qq(ibnd+i-1,j)+dqq*dxm(ibnd+i-1)
        enddo
        enddo
      endif
      
      return
      end
