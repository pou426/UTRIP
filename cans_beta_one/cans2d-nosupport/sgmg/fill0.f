c     ---------------------------------------------------------------
      subroutine fill0(u,nx)
      integer nx
      double precision u(-1:nx,-1:nx)
      integer i,j
      do j=0,nx
         do i=0,nx
            u(i,j)=0.d0
         enddo
      enddo
      return
      end
