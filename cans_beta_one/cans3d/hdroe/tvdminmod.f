c----------------------------------------------------------------------|
c     version 1.1 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|
      subroutine tvdminmod(mdir,da,daw,ix,jx,kx)

      implicit double precision (a-h,o-z)

      dimension da(ix,jx,kx)
      dimension daw(ix,jx,kx,2)

c----------------------------------------------------------------------|
c     define limiter functions
c     1. minmod limiter

      flmt(a,b)=max(0.0d0,min(b*sign(1.0d0,a),abs(a)))*sign(1.0d0,a)
c----------------------------------------------------------------------|

      if (mdir .eq. 1) then
        do k=2,kx-2
        do j=2,jx-2
        do i=2,ix-2
           daw(i,j,k,1)=da(i,j,k)
     &       +0.5d0*flmt(da(i,j,k)-da(i-1,j,k),da(i+1,j,k)-da(i,j,k))
           daw(i,j,k,2)=da(i+1,j,k)
     &       -0.5d0*flmt(da(i+2,j,k)-da(i+1,j,k),da(i+1,j,k)-da(i,j,k))
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 2) then
        do i=2,ix-2
        do k=2,kx-2
        do j=2,jx-2
           daw(i,j,k,1)=da(i,j,k)  
     &       +0.5d0*flmt(da(i,j,k)-da(i,j-1,k),da(i,j+1,k)-da(i,j,k))
           daw(i,j,k,2)=da(i,j+1,k)
     &       -0.5d0*flmt(da(i,j+2,k)-da(i,j+1,k),da(i,j+1,k)-da(i,j,k))
        enddo
        enddo
        enddo
      endif

      if (mdir .eq. 3) then
        do j=2,jx-2
        do i=2,ix-2
        do k=2,kx-2
           daw(i,j,k,1)=da(i,j,k)  
     &       +0.5d0*flmt(da(i,j,k)-da(i,j,k-1),da(i,j,k+1)-da(i,j,k))
           daw(i,j,k,2)=da(i,j,k+1)
     &       -0.5d0*flmt(da(i,j,k+2)-da(i,j,k+1),da(i,j,k+1)-da(i,j,k))
        enddo
        enddo
        enddo
      endif

      return
      end







