c======================================================================|
      subroutine bdperx(margin0,margin1,qq,ix,jx)
c======================================================================|
c
c NAME  bdperx
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

        ibnd0=1+margin0
        ibnd1=ix-margin1

        do i=1,margin0
        do j=1,jx
          qq(ibnd0-i,j) = qq(ibnd1+1-i,j)
        enddo
        enddo

        do i=1,margin1
        do j=1,jx
          qq(ibnd1+i,j) = qq(ibnd0-1+i,j)
        enddo
        enddo

      
      return
      end
