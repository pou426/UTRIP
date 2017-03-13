c======================================================================|
      subroutine pertub(vx,x,ix)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix)
      dimension vx(ix)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=1.d-6
      rkx=2.d0*2.d0*pi

      do i=1,ix
         vx(i) = vx(i)+amp*sin(rkx*x(i))
      enddo
      
      return
      end
