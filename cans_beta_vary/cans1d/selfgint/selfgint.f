c======================================================================|
      subroutine selfgint(gx,gxm,g0,ro,xmc,x,dx,dxm,ix)
c======================================================================|
c
c NAME  selfg
c
c PURPOSE
c    derive self-gravity by simply integrate mass from center of mass
c
c OUTPUTS
c    gx(ix), gxm(ix) : [double] gravitational acceleration
c
c INPUTS
c    g0: [double] strength of gravity
c    ro(ix): [double] density
c    xmc : [double] coordinate of center of mass
c    x(ix): [double] coordinate
c    dx(ix),dxm(ix): [double] grid spacing
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension gx(ix), gxm(ix)
      dimension x(ix)
      dimension dx(ix), dxm(ix)

      dimension dgx(ix)
c----------------------------------------------------------------------|      

      pi = acos(-1.0d0)
      p4=pi*4.

      do i=1,ix-1
        rk0=-dxm(i)*p4*g0*ro(i)
        rk1=-0.5*p4*g0*(ro(i)*dx(i+1)+ro(i+1)*dx(i))
        rk2=rk1
        rk3=-dxm(i)*p4*g0*ro(i+1)
        dgx(i)=1./6.*(rk0+2.*rk1+2.*rk2+rk3)
      enddo

      gx(1)=0.
      do i=2,ix
        gx(i)=gx(i-1)+dgx(i-1)
      enddo

      do i=1,ix-1
        if (xmc.ge.x(i).and.xmc.lt.x(i+1)) goto 100
      enddo
100   continue
      imc=i

      dx0=xmc-x(imc)

      gxmc=gx(imc)+dx0/dxm(imc)*dgx(imc)
      do i=1,ix
        gx(i)=gx(i)-gxmc
      enddo


      do i=1,ix-1
         gxm(i) = 0.5*(gx(i)+gx(i+1))
      enddo



      return
      end
