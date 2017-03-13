c======================================================================|
      subroutine pertub(thini,phini,vx,vy,vz,x,ix,y,jx,z,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx),z(kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=1.d-6
      rk=2.0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ss=x(i)*sin(thini)*cos(phini)+y(j)*sin(thini)*sin(phini)
     &    +z(k)*cos(thini)
         vs = amp*sin(rk*2*pi*ss)
         vx(i,j,k)  = vx(i,j,k)+vs*sin(thini)*cos(phini)
         vy(i,j,k)  = vy(i,j,k)+vs*sin(thini)*sin(phini)
         vz(i,j,k)  = vz(i,j,k)+vs*cos(thini)
      enddo
      enddo
      enddo
      
      return
      end
