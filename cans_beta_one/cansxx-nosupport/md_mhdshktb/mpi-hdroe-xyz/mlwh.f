c======================================================================|
      subroutine mlwh(uh,u2,dt
     &    ,fx,dxi,dxim,fy,dyi,dyim,fz,dzi,dzim,ix,jx,kx,mfdim)
c======================================================================|
c
c NAME  mlwhalf
c
c PURPOSE
c    first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    u2(ix,jx,kx): [double] variation in this step
c
c OUTPUTS
c    uh(ix,jx,kx) : [double] half step results on mid-grid points
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
      dimension mfdim(3)
      dimension dxi(ix),dxim(ix)
      dimension dyi(jx),dyim(jx)
      dimension dzi(kx),dzim(kx)
      dimension uh(ix,jx,kx),u2(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)
      dimension f(ix,jx,kx)
c----------------------------------------------------------------------|
c     include contribution to u2 from this step value      
c----------------------------------------------------------------------|
      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=2,ix-1
         u2(i,j,k)=u2(i,j,k)
     &             -0.5d0*dt*(0.5d0*dxi(i)*(fx(i+1,j,k)-fx(i-1,j,k)))
      enddo
      enddo
      enddo
      endif
      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=2,jx-1
      do i=1,ix
         u2(i,j,k)=u2(i,j,k)
     &             -0.5d0*dt*(0.5d0*dyi(j)*(fy(i,j+1,k)-fy(i,j-1,k)))
      enddo
      enddo
      enddo
      endif
      if (mfdim(3).eq.1) then
      do k=2,kx-1
      do j=1,jx
      do i=1,ix
         u2(i,j,k)=u2(i,j,k)
     &             -0.5d0*dt*(0.5d0*dzi(k)*(fz(i,j,k+1)-fz(i,j,k-1)))
      enddo
      enddo
      enddo
      endif
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
c  x-direction
      if (mfdim(1).eq.1) then
      if (mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        f(i,j,k)= (fx(i,j,k)+fx(i,j+1,k))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        f(i,j,k)= (fx(i,j,k)+fx(i,j,k+1))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix
        f(i,j,k)= (fx(i,j,k)+fx(i,j+1,k)
     &            +fx(i,j,k+1)+fx(i,j+1,k+1))/4.d0
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
      do i=1,ix-1
         uh(i,j,k)= uh(i,j,k)
     &     -(f(i+1,j,k)-f(i,j,k))*dxim(i)*dt
      enddo
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c  y-direction 

      if (mfdim(2).eq.1) then
      if (mfdim(3).eq.1.and.mfdim(1).eq.0) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        f(i,j,k)= (fy(i,j,k)+fy(i,j,k+1))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(3).eq.0.and.mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        f(i,j,k)= (fy(i,j,k)+fy(i+1,j,k))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(3).eq.1.and.mfdim(1).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix-1
        f(i,j,k)= (fy(i,j,k)+fy(i,j,k+1)
     &            +fy(i+1,j,k)+fy(i+1,j,k+1))/4.d0
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
      do j=1,jx-1
      do i=1,ix
         uh(i,j,k)= uh(i,j,k)
     &     -(f(i,j+1,k)-f(i,j,k))*dyim(j)*dt
      enddo
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c  z-direction 
      if (mfdim(3).eq.1) then
      if (mfdim(1).eq.1.and.mfdim(2).eq.0) then
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        f(i,j,k)= (fz(i,j,k)+fz(i+1,j,k))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(1).eq.0.and.mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        f(i,j,k)= (fz(i,j,k)+fz(i,j+1,k))/2.d0
      enddo
      enddo
      enddo
      else if (mfdim(1).eq.1.and.mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix-1
        f(i,j,k)= (fz(i,j,k)+fz(i+1,j,k)
     &            +fz(i,j+1,k)+fz(i+1,j+1,k))/4.d0
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
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
         uh(i,j,k)= uh(i,j,k)
     &     -(f(i,j,k+1)-f(i,j,k))*dzim(k)*dt
      enddo
      enddo
      enddo
      endif

      return
      end
