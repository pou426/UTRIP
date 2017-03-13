c=====================================================================|
      subroutine getfbx(fx,fy,fz,ey,ez,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension ey(ix,jx,kx),ez(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  0.d0
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fy(i,j,k)=  ez(i,j,k)
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fz(i,j,k)= -ey(i,j,k)
      enddo
      enddo
      enddo
      endif

      return
      end
