c======================================================================|
      subroutine mlwhalf(u,un,du,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
c======================================================================|
c
c NAME  mlwhalf
c
c PURPOSE
c    first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c
c OUTPUTS
c    un(ix,jx) : [double] half step results on mid-grid points
c
c INPUTS
c    u(ix,jx) : [double] basic variables    
c    f(ix,jx) : [double] flux in x-direction
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time 
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension dxi(ix),dxim(ix)
      dimension dyi(jx),dyim(jx)
      dimension u(ix,jx), un(ix,jx), du(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         du(i,j)=du(i,j)-0.5*dt*(0.5*dxi(i)*(fx(i+1,j)-fx(i-1,j)))
     &                  -0.5*dt*(0.5*dyi(j)*(fy(i,j+1)-fy(i,j-1)))
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         uh   = 0.25*(u(i+1,j)+u(i,j)+u(i+1,j+1)+u(i,j+1))
         dfdx = dxim(i)*(fx(i+1,j)-fx(i,j)+fx(i+1,j+1)-fx(i,j+1))/2
         dfdy = dyim(j)*(fy(i,j+1)-fy(i,j)+fy(i+1,j+1)-fy(i+1,j))/2
         un(i,j)= uh-dt*dfdx-dt*dfdy
      enddo
      enddo

      return
      end
