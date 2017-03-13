c======================================================================|
      subroutine mlwf(u2,dt
     &    ,fx,dxi,ux0,ux1,fy,dyi,uy0,uy1,fz,dzi,uz0,uz1,ix,jx,kx,mfdim)
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
      dimension mfdim(3)
      dimension dxi(ix),ux0(ix),ux1(ix)
      dimension dyi(jx),uy0(jx),uy1(jx)
      dimension dzi(kx),uz0(kx),uz1(kx)
      dimension u2(ix,jx,kx),fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
      dimension f(ix,jx,kx)
c----------------------------------------------------------------------|      

      if (mfdim(1).eq.1) then
      if (mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=2,jx
      do i=1,ix
        f(i,j,k)= uy1(j)*fx(i,j,k)+uy0(j)*fx(i,j-1,k)
      enddo
      enddo
      enddo
      else if (mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=2,kx
      do j=1,jx
      do i=1,ix
        f(i,j,k)= uz1(k)*fx(i,j,k)+uz0(k)*fx(i,j,k-1)
      enddo
      enddo
      enddo
      else if (mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=2,kx
      do j=2,jx
      do i=1,ix
        f(i,j,k)= uy1(j)*uz1(k)*fx(i,j,k)
     &           +uy0(j)*uz1(k)*fx(i,j-1,k)
     &           +uy1(j)*uz0(k)*fx(i,j,k-1)
     &           +uy0(j)*uz0(k)*fx(i,j-1,k-1)
      enddo
      enddo
      enddo
      else
      do k=1,kx
      do j=1,jx
      do i=1,ix
        f(i,j,k)= fx(i,j,k)
      enddo
      enddo
      enddo
      endif

      do k=1,kx
      do j=1,jx
      do i=2,ix
         u2(i,j,k)= u2(i,j,k)-0.5d0*dt*dxi(i)*(f(i,j,k)-f(i-1,j,k))
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      if (mfdim(3).eq.1.and.mfdim(1).eq.0) then
      do k=2,kx
      do j=1,jx
      do i=1,ix
        f(i,j,k)= uz1(k)*fy(i,j,k)+uz0(k)*fy(i,j,k-1)
      enddo
      enddo
      enddo
      else if (mfdim(3).eq.0.and.mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=2,ix
        f(i,j,k)= ux1(i)*fy(i,j,k)+ux0(i)*fy(i-1,j,k)
      enddo
      enddo
      enddo
      else if (mfdim(3).eq.1.and.mfdim(1).eq.1) then
      do k=2,kx
      do j=1,jx
      do i=2,ix
        f(i,j,k)= uz1(k)*ux1(i)*fy(i,j,k)
     &           +uz0(k)*ux1(i)*fy(i,j,k-1)
     &           +uz1(k)*ux0(i)*fy(i-1,j,k)
     &           +uz0(k)*ux0(i)*fy(i-1,j,k-1)
      enddo
      enddo
      enddo
      else
      do k=1,kx
      do j=1,jx
      do i=1,ix
        f(i,j,k)= fy(i,j,k)
      enddo
      enddo
      enddo
      endif

      do k=1,kx
      do j=2,jx
      do i=1,ix
         u2(i,j,k)= u2(i,j,k)-0.5d0*dt*dyi(j)*(f(i,j,k)-f(i,j-1,k))
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      if (mfdim(1).eq.1.and.mfdim(2).eq.0) then
      do k=1,kx
      do j=1,jx
      do i=2,ix
        f(i,j,k)= ux1(i)*fz(i,j,k)+ux0(i)*fz(i-1,j,k)
      enddo
      enddo
      enddo
      else if (mfdim(1).eq.0.and.mfdim(2).eq.1) then
      do k=1,kx
      do j=2,jx
      do i=1,ix
        f(i,j,k)= uy1(j)*fz(i,j,k)+uy0(j)*fz(i,j-1,k)
      enddo
      enddo
      enddo
      else if (mfdim(1).eq.1.and.mfdim(2).eq.1) then
      do k=1,kx
      do j=2,jx
      do i=2,ix
        f(i,j,k)= ux1(i)*uy1(j)*fz(i,j,k)
     &           +ux0(i)*uy1(j)*fz(i-1,j,k)
     &           +ux1(i)*uy0(j)*fz(i,j-1,k)
     &           +ux0(i)*uy0(j)*fz(i-1,j-1,k)
      enddo
      enddo
      enddo
      else
      do k=1,kx
      do j=1,jx
      do i=1,ix
        f(i,j,k)= fz(i,j,k)
      enddo
      enddo
      enddo
      endif

      do k=2,kx
      do j=1,jx
      do i=1,ix
         u2(i,j,k)= u2(i,j,k)-0.5d0*dt*dzi(k)*(f(i,j,k)-f(i,j,k-1))
      enddo
      enddo
      enddo
      endif

      return
      end
