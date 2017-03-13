c======================================================================|
      subroutine getfee(fx,fy,fz,ro,pr,eh,vx,vy,vz,bx,by,bz,ex,ey,ez
     &      ,mu,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),eh(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      double precision mu

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ep= eh(i,j,k) + pr(i,j,k)
     &     + (vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)*ro(i,j,k)/2.d0

         fx(i,j,k)= ep*vx(i,j,k)
     &              +(bz(i,j,k)*ey(i,j,k)-by(i,j,k)*ez(i,j,k))/mu
      enddo
      enddo
      enddo
      endif
      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ep= eh(i,j,k) + pr(i,j,k)
     &     + (vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)*ro(i,j,k)/2.d0

         fy(i,j,k)= ep*vy(i,j,k)
     &              +(bx(i,j,k)*ez(i,j,k)-bz(i,j,k)*ex(i,j,k))/mu
      enddo
      enddo
      enddo
      endif
      if (mfdim(3).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         ep= eh(i,j,k) + pr(i,j,k)
     &     + (vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)*ro(i,j,k)/2.d0

         fz(i,j,k)= ep*vz(i,j,k)
     &              +(by(i,j,k)*ex(i,j,k)-bx(i,j,k)*ey(i,j,k))/mu
      enddo
      enddo
      enddo
      endif

      return
      end
