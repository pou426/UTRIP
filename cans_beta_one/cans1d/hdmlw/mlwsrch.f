c======================================================================|
      subroutine mlwsrch(un,du,dt,s,ix)
c======================================================================|
c
c NAME  mlwsrch
c
c PURPOSE
c    apply source term for first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix): [double] variation in this step  
c    un(ix) : [double] half step results on mid-grid points
c
c OUTPUTS
c    None
c
c INPUTS
c    s(ix) : [double] source term
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension un(ix), du(ix)
      dimension s(ix)
c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do i=2,ix-1
         du(i)=du(i)+0.5*dt*s(i)
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do i=1,ix-1
c     cell average
         sh   = 0.5*(s(i+1)+s(i))
c     summation of all terms 
         un(i)= un(i)+dt*sh
      enddo

      return
      end
