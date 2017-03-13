c======================================================================|
      subroutine ccipadvrd(rdm,ro,dt,vx,dxm,ix)
c======================================================================|
c
c NAME  ccipadvrd
c
c PURPOSE
c    solve eqs. by Conservative CIP
c        * simple advection
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    rodx(ix): [double] density gradient
c    rdm(ix) : [double] integrated mass in a cell
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vx(ix), vxm(ix) : [double] velocity along the x-cordinate
c    dx(ix), dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dxm(ix)
      dimension vx(ix)
      dimension ro(ix)
      dimension roh(ix)
      dimension rdm(ix),rdmh(ix),drd(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix
        rdmh(i)=rdm(i)
        roh(i)=ro(i)
      enddo

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      isft=0
      do i=2,ix-1
        if (vx(i).le.0.) then
          sgnu=-1.d0
          iup=i+1
          iupm=i+isft
          dx1=+dxm(iupm)
        else
          sgnu=+1.d0
          iup=i-1
          iupm=i-1+isft
          dx1=-dxm(iupm)
        endif

        a= sgnu*6.d0*rdmh(iupm)/dx1**3 +3.d0*(roh(i)+roh(iup))/dx1**2
        b=-sgnu*6.d0*rdmh(iupm)/dx1**2 -2.d0*(2.*roh(i)+roh(iup))/dx1

        xx=-vx(i)*dt
        drd(i)=-(a/3.d0*xx**2+b/2.d0*xx+roh(i))*xx
        roh(i)=a*xx**2+b*xx+roh(i)
      enddo
      do i=1,ix-1
        rdmh(i)=rdmh(i)+drd(i)-drd(i+1)
      enddo

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do i=1,ix
        ro(i)=roh(i)
        rdm(i)=rdmh(i)
      enddo


      return
      end
