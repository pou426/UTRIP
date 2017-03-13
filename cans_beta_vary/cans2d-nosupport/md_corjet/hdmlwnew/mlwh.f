c======================================================================|
      subroutine mlwh(uh,u2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
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
c    uh(ix,jx) : [double] half step results on mid-grid points
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
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx)
      dimension uh(ix,jx), u2(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      call dx2dxi(dx,dxi,ix)
      call dx2dxi(dxm,dxim,ix)
      call dx2dxi(dy,dyi,jx)
      call dx2dxi(dym,dyim,jx)

c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         u2(i,j)=u2(i,j)
     &      -0.5d0*dt*(0.5d0*dxi(i)*(fx(i+1,j)-fx(i-1,j)))
     &      -0.5d0*dt*(0.5d0*dyi(j)*(fy(i,j+1)-fy(i,j-1)))
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         dfdx = dxim(i)*(fx(i+1,j)-fx(i,j)+fx(i+1,j+1)-fx(i,j+1))/2.d0
         dfdy = dyim(j)*(fy(i,j+1)-fy(i,j)+fy(i+1,j+1)-fy(i+1,j))/2.d0
         uh(i,j)= uh(i,j)-dt*dfdx-dt*dfdy
      enddo
      enddo

      return
      end
