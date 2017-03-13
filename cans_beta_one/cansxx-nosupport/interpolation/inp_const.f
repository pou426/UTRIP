c----------------------------------------------------------------------|
c     version 1.1 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|
      subroutine inp_const(da,daw,ix,jx,kx,mdir)

      implicit double precision (a-h,o-z)

      dimension da(ix,jx,kx)
      dimension daw(ix,jx,kx,2)

c----------------------------------------------------------------------|

      if (mdir .eq. 1) then
        do k=1,kx
        do j=1,jx
        do i=1,ix-1
           daw(i,j,k,1)=da(i,j,k)
           daw(i,j,k,2)=da(i+1,j,k)
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 2) then
        do k=1,kx
        do j=1,jx-1
        do i=1,ix
           daw(i,j,k,1)=da(i,j,k)  
           daw(i,j,k,2)=da(i,j+1,k)
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 3) then
        do k=1,kx-1
        do j=1,jx
        do i=1,ix
           daw(i,j,k,1)=da(i,j,k)  
           daw(i,j,k,2)=da(i,j,k+1)
        enddo
        enddo
        enddo
      endif

      return
      end







