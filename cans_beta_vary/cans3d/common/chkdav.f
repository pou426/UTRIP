c======================================================================|
      subroutine chkdav(n_events,da,vx,vy,vz,floor,ix,jx,kx)
c======================================================================|
c 
c NAME  chkdav
c 
c PURPOSE
c    keep values of data greater than the threshold.
c    If the value is less than the threshold,
c    set velocity to zero at the same point.
c    
c INPUTS & OUTPUTS
c    da(ix,jx,kx): [double] variable
c    vx(ix,jx,kx): [double] variable
c    vy(ix,jx,kx): [double] variable
c 
c OUTPUTS
c    n_events: [integer] # of data points which are less than threshold
c    
c INPUTS
c    floor: [double] threshold value
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension da(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
c----------------------------------------------------------------------|

      n_events=0
      
      do k=1,kx
      do j=1,jx
      do i=1,ix
         if (da(i,j,k).lt.floor) then
             n_events=n_events+1
             da(i,j,k) = floor
             vx(i,j,k) = 0.0
             vy(i,j,k) = 0.0
             vz(i,j,k) = 0.0
         endif
      enddo
      enddo
      enddo
      

      return
      end      
