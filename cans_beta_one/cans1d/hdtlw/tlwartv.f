c======================================================================|
      subroutine tlwartv(u,du,dt,qx,dxi,dxim,ix)
c======================================================================|
c
c NAME  tlwartv
c
c PURPOSE
c    apply artificial viscosity 
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step  
c
c OUTPUTS
c    None
c
c INPUTS
c    u(ix) : [double] basic variables    
c    qx(ix) : [double] coefficients of artificial viscosity 
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
      dimension u(ix),du(ix),qx(ix)
c----------------------------------------------------------------------|
      do i=2,ix-1
         du(i)=du(i)+dt*(dxi(i)*(qx(i)*dxim(i)*(u(i+1)-u(i))
     &        - qx(i-1)*dxim(i-1)*(u(i)-u(i-1))))
      enddo

      return
      end
