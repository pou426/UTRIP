c======================================================================|
      subroutine pertub(thini,vx,vy,x,ix,y,jx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx)
      dimension vx(ix,jx),vy(ix,jx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=1.d-6
      rk=2.0

      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         vs = amp*sin(rk*2*pi*ss)
         vx(i,j)  = vx(i,j)+vs*cos(thini)
         vy(i,j)  = vy(i,j)+vs*sin(thini)
      enddo
      enddo
      
      return
      end
