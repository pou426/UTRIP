c======================================================================|
      subroutine mlwhalf(u,un,du,dt,f,dxi,dxim,ix)
c======================================================================|
c
c NAME  mlwhalf
c
c PURPOSE
c    first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step
c
c OUTPUTS
c    un(ix) : [double] half step results on mid-grid points
c
c INPUTS
c    u(ix) : [double] basic variables    
c    f(ix) : [double] flux in x-direction
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxi(ix),dxim(ix)
      dimension u(ix), un(ix), du(ix), f(ix)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do i=2,ix-1
         du(i)=du(i)-0.5*dt*(0.5*dxi(i)*(f(i+1)-f(i-1)))
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do i=1,ix-1
c     cell average
         uh   = 0.5*(u(i+1)+u(i))
c     flux across cell boundary
         dfdx = dxim(i)*(f(i+1)-f(i))
c     summation of all terms 
         un(i)= uh-dt*dfdx
      enddo
c      
      return
      end
