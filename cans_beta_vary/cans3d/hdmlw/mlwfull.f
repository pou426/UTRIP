c======================================================================|
      subroutine mlwfull(du,dt
     &    ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
c======================================================================|
c
c NAME  mlwfull
c
c PURPOSE
c    second half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx,kx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    f(ix,jx,kx) : [double] flux in x-direction
c    dt: [double] delta time
c    dxi(ix), dxim(ix) : [double] 1/dx
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxi(ix),ux0(ix),ux1(ix)
      dimension dyi(jx),uy0(jx),uy1(jx)
      dimension dzi(kx),uz0(kx),uz1(kx)
      dimension du(ix,jx,kx),fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
c----------------------------------------------------------------------|      

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         dfdx = dxi(i)*
     &         (uz0(k)*(uy0(j)*(fx(i,j-1,k-1)-fx(i-1,j-1,k-1))
     &                 +uy1(j)*(fx(i,j  ,k-1)-fx(i-1,j  ,k-1)))
     &         +uz1(k)*(uy0(j)*(fx(i,j-1,k  )-fx(i-1,j-1,k  ))
     &                 +uy1(j)*(fx(i,j  ,k  )-fx(i-1,j  ,k  ))))
         dfdy = dyi(j)*
     &         (ux0(i)*(uz0(k)*(fy(i-1,j,k-1)-fy(i-1,j-1,k-1))
     &                 +uz1(k)*(fy(i-1,j,k  )-fy(i-1,j-1,k  )))
     &         +ux1(i)*(uz0(k)*(fy(i  ,j,k-1)-fy(i  ,j-1,k-1))
     &                 +uz1(k)*(fy(i  ,j,k  )-fy(i  ,j-1,k  ))))
         dfdz = dzi(k)*
     &         (uy0(j)*(ux0(i)*(fz(i-1,j-1,k)-fz(i-1,j-1,k-1))
     &                 +ux1(i)*(fz(i  ,j-1,k)-fz(i  ,j-1,k-1)))
     &         +uy1(j)*(ux0(i)*(fz(i-1,j  ,k)-fz(i-1,j  ,k-1))
     &                 +ux1(i)*(fz(i  ,j  ,k)-fz(i  ,j  ,k-1))))
         du(i,j,k)= du(i,j,k)-0.5*dt*dfdx-0.5*dt*dfdy-0.5*dt*dfdz
      enddo
      enddo
      enddo

      return
      end
