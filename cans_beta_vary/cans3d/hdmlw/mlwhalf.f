c======================================================================|
      subroutine mlwhalf(u,un,du,dt
     &    ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
c======================================================================|
c
c NAME  mlwhalf
c
c PURPOSE
c    first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx,kx): [double] variation in this step
c
c OUTPUTS
c    un(ix,jx,kx) : [double] half step results on mid-grid points
c
c INPUTS
c    u(ix,jx,kx) : [double] basic variables    
c    f(ix,jx,kx) : [double] flux in x-direction
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time 
c    ix,jx,kx: [integer] dimension size
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension dxi(ix),dxim(ix)
      dimension dyi(jx),dyim(jx)
      dimension dzi(kx),dzim(kx)
      dimension  u(ix,jx,kx),un(ix,jx,kx),du(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         du(i,j,k)=du(i,j,k)
     &             -0.5*dt*(0.5*dxi(i)*(fx(i+1,j,k)-fx(i-1,j,k)))
     &             -0.5*dt*(0.5*dyi(j)*(fy(i,j+1,k)-fy(i,j-1,k)))
     &             -0.5*dt*(0.5*dzi(k)*(fz(i,j,k+1)-fz(i,j,k-1)))
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        uh = (u(i+1,j,k  )+u(i,j,k  )+u(i+1,j+1,k  )+u(i,j+1,k  ) 
     &       +u(i+1,j,k+1)+u(i,j,k+1)+u(i+1,j+1,k+1)+u(i,j+1,k+1))/8
      dfdx = (fx(i+1,j,k  )-fx(i,j,k  )+fx(i+1,j+1,k  )-fx(i,j+1,k  )
     &       +fx(i+1,j,k+1)-fx(i,j,k+1)+fx(i+1,j+1,k+1)-fx(i,j+1,k+1))
     &       /4*dxim(i)
      dfdy = (fy(i,j+1,k  )-fy(i,j,k  )+fy(i+1,j+1,k  )-fy(i+1,j,k  )
     &       +fy(i,j+1,k+1)-fy(i,j,k+1)+fy(i+1,j+1,k+1)-fy(i+1,j,k+1))
     &       /4*dyim(j)
      dfdz = (fz(i,j  ,k+1)-fz(i,j  ,k)+fz(i+1,j  ,k+1)-fz(i+1,j  ,k)
     &       +fz(i,j+1,k+1)-fz(i,j+1,k)+fz(i+1,j+1,k+1)-fz(i+1,j+1,k))
     &       /4*dzim(k)
         un(i,j,k)= uh-dt*dfdx-dt*dfdy-dt*dfdz
      enddo
      enddo
      enddo

      return
      end
