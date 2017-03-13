c======================================================================|
      subroutine mlwsrch(un,du,dt,s,ix,jx,kx)
c======================================================================|
c
c NAME  mlwsrch
c
c PURPOSE
c    apply source term for first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx,kx): [double] variation in this step
c    un(ix,jx,kx) : [double] half step results on mid-grid points
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
      dimension un(ix,jx,kx), du(ix,jx,kx)
      dimension s(ix,jx,kx)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         du(i,j,k)=du(i,j,k)+0.5*dt*s(i,j,k)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        sh = (s(i+1,j,k  )+s(i,j,k  )+s(i+1,j+1,k  )+s(i,j+1,k  )
     &       +s(i+1,j,k+1)+s(i,j,k+1)+s(i+1,j+1,k+1)+s(i,j+1,k+1))/8
        un(i,j,k)= un(i,j,k)+dt*sh
      enddo
      enddo
      enddo

      return
      end
