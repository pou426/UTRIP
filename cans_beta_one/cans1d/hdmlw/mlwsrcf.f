c======================================================================|
      subroutine mlwsrcf(du,dt,s,ux0,ux1,ix)
c======================================================================|
c
c NAME  mlwsrcf
c
c PURPOSE
c    apply source term for second half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step  
c
c OUTPUTS
c    None
c
c INPUTS
c    s(ix) : [double] source term
c    ux0(ix) : [double] 0.5*dx(i)/dxm(i)
c    ux1(ix) : [double] 0.5*dx(i-1)/dxm(i)
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension ux0(ix),ux1(ix)
      dimension du(ix),s(ix) 
c----------------------------------------------------------------------|      
      do i=2,ix-1
         sh   = ux0(i)*s(i-1)+ux1(i)*s(i)
         du(i)= du(i)+0.5*dt*sh
      enddo

      return
      end
