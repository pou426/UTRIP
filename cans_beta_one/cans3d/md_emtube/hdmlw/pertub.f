c======================================================================|
      subroutine pertub(ro,x,ix,y,jx,z,kx,rtube,ytube,ztube)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx),z(kx)
      dimension ro(ix,jx,kx)
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      amp=0.1d0
      wptb=25.d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
        rr=sqrt(y(j)**2+(z(k)-ztube)**2)
        if (rr.le.rtube) then
        ro(i,j,k) = ro(i,j,k)*(1-
     &     amp*cos(2*pi*x(i)/wptb)
     &    *0.5*(tanh((x(i)+0.75*wptb)/0.5)-tanh((x(i)-0.75*wptb)/0.5)))
        endif
      enddo
      enddo
      enddo

      return
      end
