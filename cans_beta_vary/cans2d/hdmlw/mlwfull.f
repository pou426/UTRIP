c======================================================================|
      subroutine mlwfull(du,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
c======================================================================|
c
c NAME  mlwfull
c
c PURPOSE
c    second half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    f(ix,jx) : [double] flux in x-direction
c    dt: [double] delta time
c    dxi(ix), dxim(ix) : [double] 1/dx
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxi(ix),ux0(ix),ux1(ix)
      dimension dyi(jx),uy0(jx),uy1(jx)
      dimension du(ix,jx),fx(ix,jx),fy(ix,jx)
c----------------------------------------------------------------------|      
      do j=2,jx-1
      do i=2,ix-1
         dfdx = dxi(i)*
     &          (uy0(j)*(fx(i,j-1)-fx(i-1,j-1))
     &          +uy1(j)*(fx(i,j)  -fx(i-1,j)  ))
         dfdy = dyi(j)*
     &          (ux0(i)*(fy(i-1,j)-fy(i-1,j-1))
     &          +ux1(i)*(fy(i,j)  -fy(i,j-1)  ))
         du(i,j)= du(i,j)-0.5*dt*dfdx-0.5*dt*dfdy
      enddo
      enddo

      return
      end
