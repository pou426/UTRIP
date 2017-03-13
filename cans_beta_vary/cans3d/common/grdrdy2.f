c======================================================================|
      subroutine grdrdy2(dx,xm,dxm,x,ix)
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
      dimension dx(ix),x(ix)
      dimension dxm(ix),xm(ix)
c----------------------------------------------------------------------|

      do i=2,ix
         dx(i)  = 0.5d0*(dxm(i-1)+dxm(i))
      enddo
      dx(1)=dxm(1)

      do i=1,ix-1
         xm(i)=0.5d0*(x(i)+x(i+1))
      enddo
      xm(ix)=x(ix)+0.5d0*dxm(ix)

      return
      end
