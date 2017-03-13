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

      amp=0.1
      xptb=20.
      yptb1=-4.
      yptb2= 4.
      wptb=0.5

      do j=1,jx
      do i=1,ix
        dycen=abs(y(j))
        vy(i,j) = vy(i,j)
     &    +amp*cos(2*pi*x(i)/xptb)
     &   *0.5*(tanh((x(i)+0.75*xptb)/wptb)-tanh((x(i)-0.75*xptb)/wptb))
     &    *0.5*(tanh((y(j)-yptb1)/wptb)-tanh((y(j)-yptb2)/wptb))
      enddo
      enddo


      return
      end









