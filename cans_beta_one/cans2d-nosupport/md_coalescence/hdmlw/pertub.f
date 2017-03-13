c======================================================================|
      subroutine pertub(vx,x,ix,y,jx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix), y(jx)
      dimension vx(ix,jx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=0.05
      yptb1=-1.
      yptb2=-yptb1

      do j=1,jx
      do i=1,ix
        vx(i,j) = vx(i,j)
     &    -amp*sin(x(i)/2)
     &    *0.5*(tanh((y(j)-yptb1)/0.5)-tanh((y(j)-yptb2)/0.5))
      enddo
      enddo

      return
      end
