c======================================================================|
      subroutine mlwsrch1(un,du,dt,s,ro,gym,ix,jx)
cNOT YET COMPLETED! DONOT USE!
c======================================================================|
c
c NAME  mlwsrch1
c
c PURPOSE
c    apply source term for first half of Modified Lax-Wendroff method
c    for convection problem
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
c    modified 2005-7-4 H. Isobe to match the BC for convection problem
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension un(ix,jx), du(ix,jx)
      dimension s(ix,jx)
      dimension ro(ix,jx), gym(ix,jx)

c----------------------------------------------------------------------|
c     include contribution to du from this step value      
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         du(i,j)=du(i,j)+0.5*dt*s(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step using flux across cell boundary  
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
c     cell average
         rh   = (ro(i+1,j)+ro(i,j)+ro(i+1,j+1)+ro(i,j+1))/4*gym(i,j)
c     summation of all terms 
c         un(i,j)= un(i,j)+dt*sh
         un(i,j)= un(i,j)+dt*rh

c for debug
c         if (i .eq. 10) then 
c            if (j .eq. 3) then
c               write(6,*) 'srch1',gym(i,j),dt,rh
c            endif
c         endif
      enddo
      enddo

      return
      end
