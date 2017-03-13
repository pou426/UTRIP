c======================================================================|
      subroutine getfeb(fx,fy,fz,bx,by,bz,ex,ey,ez
     &      ,mu,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      double precision mu

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= +(bz(i,j,k)*ey(i,j,k)-by(i,j,k)*ez(i,j,k))/mu
      enddo
      enddo
      enddo
      endif
      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fy(i,j,k)= +(bx(i,j,k)*ez(i,j,k)-bz(i,j,k)*ex(i,j,k))/mu
      enddo
      enddo
      enddo
      endif
      if (mfdim(3).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fz(i,j,k)= +(by(i,j,k)*ex(i,j,k)-bx(i,j,k)*ey(i,j,k))/mu
      enddo
      enddo
      enddo
      endif

      return
      end
