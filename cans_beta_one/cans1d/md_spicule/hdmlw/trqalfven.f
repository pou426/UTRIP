c***********************************************************************
      subroutine trqalfven(vy,dt,time,x,ix,rr)
c***********************************************************************
c
c  torque input
c
c***********************************************************************
      implicit double precision (a-h,o-z)

      dimension x(ix)
      dimension vy(ix)
      dimension trq(ix)
      dimension rr(ix)
      dimension q(ix)
c----------------------------------------------------------------------|
c  torque condition
c----------------------------------------------------------------------|
c  dztrq : distance from loop top to flare site
      pi = acos(-1.0d0)

c----------------------------------------------------------------------|

      do i=1,ix 
         call random_number(anum)
         q(i)=2*rr(i)*(anum-0.5)
         trq(i)=tanh((x(i)-0.75)/0.075)-1
      enddo
         trq(ix)=trq(ix-2)
         q(ix)=q(ix-2)
c======================================================================@
c  wave torque input
c======================================================================@
      trq0=30.


      do i=1,ix
        vy(i)=vy(i)+dt*trq0*q(i)*trq(i)
      enddo


      return
      end
