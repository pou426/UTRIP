c======================================================================|
      subroutine bdpery(margin0,margin1,qq,ix,jx)
c======================================================================|
c
c NAME  bdpery
c
c PURPOSE
c    apply periodic boundary condition
c
c INPUTS & OUTPUTS
c    qq(ix,jx): [double] variable
c
c OUTPUTS
c    None
c
c INPUTS
c    ix,jx: [integer] dimension size
c    margin0: [integer] margin, i.e. # of grid points outside the boundary
c                   for smaller i side.
c    margin1: [integer] margin, i.e. # of grid points outside the boundary
c                   for larger i side.
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension qq(ix,jx)
c----------------------------------------------------------------------|

        jbnd0=1+margin0
        jbnd1=jx-margin1

        do j=1,margin0
        do i=1,ix
          qq(i,jbnd0-j) = qq(i,jbnd1+1-j)
        enddo
        enddo

        do j=1,margin1
        do i=1,ix
          qq(i,jbnd1+j) = qq(i,jbnd0-1+j)
        enddo
        enddo

      
      return
      end
