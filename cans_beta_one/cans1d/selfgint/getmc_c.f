c======================================================================|
      subroutine getmc_c(xmc,ro,dv,x,dx,ix)
c======================================================================|
c
c NAME  getmc_c
c
c PURPOSE
c    get center of mass
c       * non-uniform corss section
c
c OUTPUTS
c    xmc : [double] coordinate of center of mass
c
c INPUTS
c    ro(ix): [double] density
c    x(ix): [double] coordinate
c    dx(ix): [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension ro(ix)
      dimension x(ix),dx(ix)
      dimension dv(ix)
c----------------------------------------------------------------------|      


      rmom=0.
      rmas=0.
      do i=1,ix
        rmas=rmas+dv(i)*ro(i)*dx(i)
        rmom=rmom+dv(i)*ro(i)*dx(i)*x(i)
      enddo

      xmc=rmom/rmas

      return
      end
