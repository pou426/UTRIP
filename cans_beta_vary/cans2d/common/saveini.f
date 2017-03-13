c======================================================================|
      subroutine saveini(daini,da,ix,jx)
c======================================================================|
c
c NAME  saveini
c
c PURPOSE
c    save initial state for boundary conditions etc.
c
c OUTPUTS
c    daini(ix,jx): [double] variable in initial condition
c
c INPUTS
c    da(ix,jx): [double] variable
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension da(ix,jx),daini(ix,jx)
c----------------------------------------------------------------------|

      do i=1,ix
      do j=1,jx
        daini(i,j)=da(i,j)
      enddo
      enddo

      return
      end
