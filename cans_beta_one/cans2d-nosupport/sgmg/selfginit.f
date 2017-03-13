c***** self-gravity only first step
      subroutine selfginit(uod,nx)
      double precision uod(-1:nx,-1:nx)
      do j=-1,nx
         do i=-1,nx
            uod(i,j)=0.0d0
         enddo
      enddo
      return
      end
