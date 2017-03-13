c======================================================================|
      subroutine getfry(fx,fy,ro,pr,vx,vy,bx,by,mu,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      double precision mu

      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)-by(i,j)*bx(i,j)/mu
         fy(i,j)= ro(i,j)*vy(i,j)**2+pr(i,j)
     &          +(bx(i,j)**2-by(i,j)**2)/(2.d0*mu)
      enddo
      enddo

      return
      end
