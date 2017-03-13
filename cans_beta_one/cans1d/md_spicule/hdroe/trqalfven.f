c***********************************************************************
      subroutine trqalfven(vy,dt,time,x,ix)
c***********************************************************************
c
c  torque input
c
c***********************************************************************
      implicit double precision (a-h,o-z)

      dimension x(ix)
      dimension vy(ix)
      dimension trq(ix)
c----------------------------------------------------------------------|
c  torque condition
c----------------------------------------------------------------------|
c  dztrq : distance from loop top to flare site
      pi = acos(-1.0d0)

      wtrq=0.75
      ztrq=0.d0
c----------------------------------------------------------------------|

      p2=2*pi

      do i=1,ix
         trq(i)=exp(-0.5*((x(i)-ztrq)/wtrq)**2)
      enddo
         trq(ix)=trq(ix-2)
c======================================================================@
c  wave torque input
c======================================================================@
      ptrq=10.
      trq0=1.

      q=sin(p2*time/ptrq)

      do i=1,ix
        vy(i)=vy(i)+dt*trq0*q*trq(i)
      enddo


      return
      end
