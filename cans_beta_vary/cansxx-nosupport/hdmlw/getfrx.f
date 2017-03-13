c======================================================================|
      subroutine getfrx(fx,fy,fz,ro,pr,vx,vy,vz,bx,by,bz,mu
     &    ,ix,jx,kx,mfdim)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension mfdim(3)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension fx(ix,jx,kx),fy(ix,jx,kx),fz(ix,jx,kx)

      double precision mu

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
     &          +(by(i,j,k)**2+bz(i,j,k)**2-bx(i,j,k)**2)/(2.d0*mu)
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fy(i,j,k)= ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
     &                       -bx(i,j,k)*by(i,j,k)/mu
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fz(i,j,k)= ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
     &                       -bx(i,j,k)*bz(i,j,k)/mu
      enddo
      enddo
      enddo
      endif

      return
      end
