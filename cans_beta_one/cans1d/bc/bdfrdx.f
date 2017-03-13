c======================================================================|
      subroutine bdfrdx(mbnd,margin,qq,dxm,ix)
c======================================================================|
c
c NAME  bdfrdx
c
c PURPOSE
c    apply free boundary condition
c    values are extended to have constant gradient
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
      dimension qq(ix),dxm(ix)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        dqq=(qq(ibnd+1)-qq(ibnd))/dxm(ibnd)
        do i=1,margin 
          qq(ibnd-i) = qq(ibnd-i+1)-dqq*dxm(ibnd-i)
        enddo
      else
        ibnd=ix-margin
        dqq=(qq(ibnd)-qq(ibnd-1))/dxm(ibnd-1)
        do i=1,margin
          qq(ibnd+i) = qq(ibnd+i-1)+dqq*dxm(ibnd+i-1)
        enddo
      endif
      
      return
      end
