c======================================================================|
      subroutine mlwsh(uh,u2,dt,s,ix,jx,kx,mfdim)
c======================================================================|
c
c NAME  mlwsrch
c
c PURPOSE
c    apply source term for first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    u2(ix,jx,kx): [double] variation in this step
c    uh(ix,jx,kx) : [double] half step results on mid-grid points
c 
c OUTPUTS
c    None
c 
c INPUTS
c    s(ix,jx,kx) : [double] source term
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension uh(ix,jx,kx), u2(ix,jx,kx)
      dimension s(ix,jx,kx)
c----------------------------------------------------------------------|
c     include contribution to u2 from this step value      
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         u2(i,j,k)=u2(i,j,k)+0.5d0*dt*s(i,j,k)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell bouhdary  
c----------------------------------------------------------------------|
c     do k=1,kx-1
c     do j=1,jx-1
c     do i=1,ix-1
c       sh = (s(i+1,j,k  )+s(i,j,k  )+s(i+1,j+1,k  )+s(i,j+1,k  )
c    &       +s(i+1,j,k+1)+s(i,j,k+1)+s(i+1,j+1,k+1)+s(i,j+1,k+1))/8.d0
c       uh(i,j,k)= uh(i,j,k)+dt*sh
c     enddo
c     enddo
c     enddo

c----------------------------------------------------------------------|
c  3D-xyz

      if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        uh(i,j,k) =uh(i,j,k)
     &   +dt*(s(i+1,j,k  )+s(i,j,k  )+s(i+1,j+1,k  )+s(i,j+1,k  )
     &       +s(i+1,j,k+1)+s(i,j,k+1)+s(i+1,j+1,k+1)+s(i,j+1,k+1))/8.d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c  2D-xy

      else if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix-1
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i+1,j,k)
     &           +  s(i,j+1,k)+s(i+1,j+1,k))/4.d0
      enddo
      enddo
      enddo

c  2D-zx

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix-1
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i+1,j,k)
     &           +  s(i,j,k+1)+s(i+1,j,k+1))/4.d0
      enddo
      enddo
      enddo

c  2D-yz

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i,j+1,k)
     &           +  s(i,j,k+1)+s(i,j+1,k+1))/4.d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c  1D-x

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i+1,j,k))/2.d0
      enddo
      enddo
      enddo

c  1D-y

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i,j+1,k))/2.d0
      enddo
      enddo
      enddo

c  1D-z

      else if (mfdim(1).eq.0.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        uh(i,j,k)= uh(i,j,k)
     &   +dt*(s(i,j,k)+s(i,j,k+1))/2.d0
      enddo
      enddo
      enddo

      endif

      return
      end
