c======================================================================|
      subroutine pertub(vz,x,ix,y,jx,z,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx),z(kx)
      dimension vz(ix,jx,kx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=0.05
      xptb=20.
      zptb1=-2.
      zptb2= 0.
      wptb=0.5

      do k=1,kx
      do j=1,jx
      do i=1,ix
        vz(i,j,k) = vz(i,j,k)
     &    +amp*cos(2*pi*x(i)/xptb)
     &    *0.5*(tanh((x(i)+0.75*xptb)/wptb)-tanh((x(i)-0.75*xptb)/wptb))
     &    *0.5*(tanh((z(k)-zptb1)/wptb)-tanh((z(k)-zptb2)/wptb))
      enddo
      enddo
      enddo

      return
      end
