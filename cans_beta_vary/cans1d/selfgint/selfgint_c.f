c======================================================================|
      subroutine selfgint_c(gx,gxm,g0,ro,xmc,dv,dvm
     &               ,x,dx,dxm,ix)
c======================================================================|
c
c NAME  selfg_c
c
c PURPOSE
c    derive self-gravity by simply integrate mass from center of mass
c       * non-uniform corss section
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
      dimension gx(ix),gxm(ix)
      dimension x(ix)
      dimension dx(ix),dxm(ix)
      dimension dgxcs(ix)
      dimension gxcsm(ix)
      dimension dv(ix),dvm(ix)
c----------------------------------------------------------------------|      

      pi = acos(-1.0d0)
      p4=pi*4.

      do i=1,ix
        dgxcs(i)=-dx(i)*p4*g0*ro(i)*dv(i)
      enddo

      gxcsm(1)=0.
      do i=2,ix
        gxcsm(i)=gxcsm(i-1)+dgxcs(i-1)
      enddo

      xmc=0.d0

       do i=1,ix-1
         if (xmc.ge.x(i).and.xmc.lt.x(i+1)) goto 100
       enddo
  100   continue
       imc=i
 
       dxmc0=xmc-x(imc)
 
       gxcsmc=gxcsm(imc)+dxmc0/dxm(imc)*dgxcs(imc)
       do i=1,ix
         gxcsm(i)=gxcsm(i)-gxcsmc
       enddo

      do i=1,ix
        gxm(i)=gxcsm(i)/dvm(i)
      enddo

      do i=1,ix-1
         gx(i) = (gxm(i)+gxm(i+1))/2
      enddo


      return
      end
