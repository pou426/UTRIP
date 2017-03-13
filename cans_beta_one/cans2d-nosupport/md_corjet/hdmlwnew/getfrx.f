c======================================================================|
      subroutine getfrx(fx,fy,ro,pr,vx,vy,bx,by,mu,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      double precision mu

      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)
     &          +(by(i,j)**2-bx(i,j)**2)/(2.d0*mu)
         fy(i,j)= ro(i,j)*vx(i,j)*vy(i,j)-bx(i,j)*by(i,j)/mu
      enddo
      enddo

      return
      end
