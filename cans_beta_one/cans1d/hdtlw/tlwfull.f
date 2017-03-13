c======================================================================|
      subroutine tlwfull(du,dt,f,dxi,ix)
c======================================================================|
c
c NAME  tlwfull
c
c PURPOSE
c    second half of Two-step Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step  
c
c OUTPUTS
c    None
c
c INPUTS
c    f(ix) : [double] flux in x-direction
c    dt: [double] delta time
c    dxi(ix), dxim(ix) : [double] 1/dx
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxi(ix)
      dimension du(ix),f(ix)
c----------------------------------------------------------------------|      
      do i=2,ix-1
         du(i)= du(i)-dt*dxi(i)*(f(i)-f(i-1))
      enddo

      return
      end
