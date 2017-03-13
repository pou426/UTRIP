c======================================================================|
      subroutine getfee(fx,fy,ro,pr,vx,vy,bx,by,ez,gm,mu,ix,jx)
c======================================================================|

      implicit double precision (a-h,o-z)

      dimension ro(ix,jx),pr(ix,jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension ez(ix,jx)
      dimension fx(ix,jx),fy(ix,jx)

      double precision mu

      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vy(i,j)**2
         eth= pr(i,j)/(gm-1.d0) 
         ep = eth + pr(i,j)+0.5d0*ro(i,j)*vv
         fx(i,j)= ep*vx(i,j) - by(i,j)*ez(i,j)/mu
         fy(i,j)= ep*vy(i,j) + bx(i,j)*ez(i,j)/mu
      enddo
      enddo

      return
      end
