c======================================================================|
      subroutine fx2c(fy,ss,fx,x,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension x(ix)
      dimension fx(ix,jx,kx),fy(ix,jx,kx)
      dimension ss(ix,jx,kx)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss(i,j,k)= -fx(i,j,k)/x(i)
      enddo
      enddo
      enddo

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fy(i,j,k)= fy(i,j,k)/x(i)
      enddo
      enddo
      enddo
      endif

      return
      end
