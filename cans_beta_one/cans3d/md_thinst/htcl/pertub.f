c======================================================================|
      subroutine pertub(pr,x,ix,y,jx,z,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix), y(jx), z(kx)
      dimension pr(ix,jx,kx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      wexp=0.3d0
      prexp=-0.1d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         r2=x(i)**2+y(j)**2+z(k)**2
         pr(i,j,k) = pr(i,j,k)*(1+prexp*exp(-r2/wexp**2))
      enddo
      enddo
      enddo
      
      return
      end
