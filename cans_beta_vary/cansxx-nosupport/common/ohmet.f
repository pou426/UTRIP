c======================================================================|
      subroutine ohmet(ex,ey,ez,bx,by,bz,vx,vy,vz,et,dx,dy,dz
     &        ,ix,jx,kx,mfdim)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension dx(ix),dy(jx),dz(kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension cx(ix,jx,kx),cy(ix,jx,kx),cz(ix,jx,kx)
      dimension et(ix,jx,kx)

c----------------------------------------------------------------------|

      call bbtocc(cx,cy,cz,bx,by,bz,dx,dy,dz,ix,jx,kx,mfdim)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ex(i,j,k) = -vy(i,j,k)*bz(i,j,k)+vz(i,j,k)*by(i,j,k)
     &              +  et(i,j,k)*cx(i,j,k)
         ey(i,j,k) = -vz(i,j,k)*bx(i,j,k)+vx(i,j,k)*bz(i,j,k)
     &              +  et(i,j,k)*cy(i,j,k)
         ez(i,j,k) = -vx(i,j,k)*by(i,j,k)+vy(i,j,k)*bx(i,j,k)
     &              +  et(i,j,k)*cz(i,j,k)
      enddo
      enddo
      enddo


      return
      end
