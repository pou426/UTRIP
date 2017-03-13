c======================================================================|
      subroutine tlwhalf(u,uh,dt,f,dxi,dxim,ix)
c======================================================================|
c
c NAME  tlwhalf
c
c PURPOSE
c    first half of Two-step Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step
c
c OUTPUTS
c    uh(ix) : [double] half step results on mid-grid points
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
      dimension u(ix), uh(ix), f(ix)
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do i=1,ix-1
         uh(i)= 0.5d0*((u(i+1)+u(i))-dt*dxim(i)*(f(i+1)-f(i)))
      enddo

      return
      end
