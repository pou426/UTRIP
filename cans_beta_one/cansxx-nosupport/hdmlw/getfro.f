c======================================================================|
      subroutine getfro(fx,fy,fz,ro,vx,vy,vz,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension ro(ix,jx,kx),vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)
      enddo
      enddo
      enddo
      endif

      return
      end
