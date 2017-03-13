c======================================================================|
      subroutine mlwreset(u2,uh,u,ix,jx,kx,mfdim)
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
      dimension mfdim(3)
      dimension u(ix,jx,kx), uh(ix,jx,kx), u2(ix,jx,kx)

      do k=1,kx
      do j=1,jx
      do i=1,ix
        u2(i,j,k)=u(i,j,k)
      enddo
      enddo
      enddo

      if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        uh(i,j,k) 
     &     = (u(i+1,j,k  )+u(i,j,k  )+u(i+1,j+1,k  )+u(i,j+1,k  )
     &       +u(i+1,j,k+1)+u(i,j,k+1)+u(i+1,j+1,k+1)+u(i,j+1,k+1))/8.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.1.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix-1
        uh(i,j,k)= (u(i,j,k)+u(i+1,j,k)
     &           +  u(i,j+1,k)+u(i+1,j+1,k))/4.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix-1
        uh(i,j,k)= (u(i,j,k)+u(i+1,j,k)
     &           +  u(i,j,k+1)+u(i+1,j,k+1))/4.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix
        uh(i,j,k)= (u(i,j,k)+u(i,j+1,k)
     &           +  u(i,j,k+1)+u(i,j+1,k+1))/4.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.1.and.mfdim(2).eq.0.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        uh(i,j,k)= (u(i,j,k)+u(i+1,j,k))/2.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.0.and.mfdim(2).eq.1.and.mfdim(3).eq.0) then
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        uh(i,j,k)= (u(i,j,k)+u(i,j+1,k))/2.d0
      enddo
      enddo
      enddo

      else if (mfdim(1).eq.0.and.mfdim(2).eq.0.and.mfdim(3).eq.1) then
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        uh(i,j,k)= (u(i,j,k)+u(i,j,k+1))/2.d0
      enddo
      enddo
      enddo

      endif


      return
      end
