c======================================================================|
      subroutine chkda(n_events,da,floor,ix,jx,kx)
c======================================================================|
c 
c NAME  chkda
c 
c PURPOSE
c    keep values of data greater than the threshold.
c        
c INPUTS & OUTPUTS
c    da(ix,jx): [double] variable
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
c----------------------------------------------------------------------|

      n_events=0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         if(da(i,j,k).lt.floor) then
           n_events=n_events+1
           da(i,j,k)=floor
         endif
      enddo
      enddo
      enddo

       return
       end
