c======================================================================|
      subroutine grdrdy(dx,xm,dxm,x,ix)
c======================================================================|
c
c NAME  grdrdy
c
c PURPOSE
c    calculate coordinate of mid-grid points and
c    grid spacing on the grid points
c
c OUTPUTS
c    dx(ix),dxm(ix): [double] grid spacing
c    xm(ix): [double] coordinate
c
c INPUTS
c    x(ix): [double] coordinate
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),xm(ix)
      dimension dxm(ix),x(ix)
c----------------------------------------------------------------------|

      do i=1,ix-1
         dxm(i)=x(i+1)-x(i)
      enddo
      dxm(ix)=dxm(ix-1)

      do i=2,ix-1
         dx(i)  = (dxm(i-1)+dxm(i))/2.d0
      enddo
      dx(1)=dx(2)
      dx(ix)=dx(ix-1)

      do i=1,ix-1
         xm(i)=(x(i)+x(i+1))/2.d0
      enddo
      xm(ix)=xm(ix-1)+dx(ix-1)


      return
      end
