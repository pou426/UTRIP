c======================================================================|
      subroutine pertub(pr,x,ix,y,jx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix), y(jx)
      dimension pr(ix,jx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      wexp=0.3d0
      prexp=-0.1d0

      do j=1,jx
      do i=1,ix
         r2=x(i)**2+y(j)**2
         pr(i,j) = pr(i,j)*(1+prexp*exp(-r2/wexp**2))
      enddo
      enddo
      
      return
      end
