c======================================================================|
      subroutine bdfrdy(mbnd,margin,qq,dym,ix,jx)
c======================================================================|
c
c NAME  bdfrdy
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
      dimension qq(ix,jx),dym(jx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        jbnd=1+margin
        do j=1,margin 
        do i=1,ix 
          dqq=(qq(i,jbnd+1)-qq(i,jbnd))/dym(jbnd)
          qq(i,jbnd-j) = qq(i,jbnd-j+1)-dqq*dym(jbnd-j)
        enddo
        enddo
      else
        jbnd=jx-margin
        do j=1,margin
        do i=1,ix 
          dqq=(qq(i,jbnd)-qq(i,jbnd-1))/dym(jbnd-1)
          qq(i,jbnd+j) = qq(i,jbnd+j-1)+dqq*dym(jbnd+j-1)
        enddo
        enddo
      endif
      
      return
      end
