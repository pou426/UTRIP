c======================================================================|
      subroutine scrdy(dsc,dscm,sc,scm,dx,dxm,ix)
c======================================================================|
c
c NAME  scrrdy
c
c PURPOSE
c    calculate cross section derivatives
c
c OUTPUTS
c    dsc(ix), dscm(ix) : [double] cross section derivative
c
c INPUTS
c    sc(ix), scm(ix) : [double] cross section
c    dx(ix),dxm(ix): [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dsc(ix),dscm(ix)
      dimension sc(ix),scm(ix)
      dimension dx(ix),dxm(ix)
c----------------------------------------------------------------------|

      do i=2,ix-1
          dsc(i) = (scm(i) - scm(i-1))/dx(i)
      enddo

      do i=2,ix-2
        dscm(i)= (sc(i+1) - sc(i))/dxm(i)
      enddo

      return
      end
