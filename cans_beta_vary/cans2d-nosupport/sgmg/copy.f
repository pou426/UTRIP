c     -------------------------------------------------------------------
      subroutine copy(aout,ain,nx)
      integer nx
      double precision ain(-1:nx,-1:nx),aout(-1:nx,-1:nx)
      integer i,j
      do j=0,nx
         do i=0,nx
            aout(i,j)=ain(i,j)
         enddo
      enddo
      return
      end
