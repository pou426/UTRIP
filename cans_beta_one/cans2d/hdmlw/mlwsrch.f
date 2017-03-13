c======================================================================|
      subroutine mlwsrch(un,du,dt,s,ix,jx)
c======================================================================|
c
c NAME  mlwsrch
c
c PURPOSE
c    apply source term for first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c    un(ix,jx) : [double] half step results on mid-grid points
c 
c OUTPUTS
c    None
c 
c INPUTS
c    s(ix,jx) : [double] source term
c    dt: [double] delta time
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension un(ix,jx), du(ix,jx)
      dimension s(ix,jx)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         du(i,j)=du(i,j)+0.5d0*dt*s(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
c     cell average
         sh   = (s(i+1,j)+s(i,j)+s(i+1,j+1)+s(i,j+1))/4.d0
c     summation of all terms 
         un(i,j)= un(i,j)+dt*sh
      enddo
      enddo

      return
      end
