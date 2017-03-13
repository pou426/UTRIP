c======================================================================|
      subroutine pertub(vy,x,ix,y,jx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix), y(jx)
      dimension vy(ix,jx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=0.05d0
      xptb=20.d0
      yptb1=-2.d0
      yptb2= 0.d0
      wptb=0.5d0

      do j=1,jx
      do i=1,ix
        vy(i,j) = vy(i,j)
     &    +amp*cos(2.d0*pi*x(i)/xptb)
     &      *0.5d0*(tanh((x(i)+0.75d0*xptb)/wptb)
     &             -tanh((x(i)-0.75d0*xptb)/wptb))
     &    *0.5d0*(tanh((y(j)-yptb1)/wptb)-tanh((y(j)-yptb2)/wptb))
      enddo
      enddo

      return
      end
