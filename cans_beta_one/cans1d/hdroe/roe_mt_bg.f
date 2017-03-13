c======================================================================|
      subroutine roe_mt_bg(ro,vx,vy,by,bx,bxm,dt,cs2
     &            ,gx,dsc,scm,dv,rr,rrm,drr,dx,ix)
c======================================================================|
c
c NAME  roe_mt_bg
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * isothermal MHD
c        * axial symmetry
c        * non-uniform poloidal magnetic field
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    by(ix): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    gx(ix), gxm(ix) : [double] gravity
c    scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    rr(ix), rrm(ix) : [double] distance from rotation axis
c    drr(ix), drrm(ix) : [double] distance gradient from rotation axis
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2005-2-4 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)
      dimension ro(ix),vx(ix),vy(ix),by(ix)
      dimension bx(ix),bxm(ix)
      dimension dsc(ix),scm(ix),dv(ix)
      dimension rosc(ix),rxsc(ix),rysc(ix),bysc(ix)
      dimension rosch(ix),rxsch(ix),rysch(ix),bysch(ix)
      dimension fro(ix),frx(ix),fry(ix),fby(ix)
      dimension roh(ix),vxh(ix),vyh(ix),byh(ix)
      dimension row(ix,2),vxw(ix,2),vyw(ix,2)
     &         ,bxw(ix,2),byw(ix,2)
      dimension gx(ix)
      dimension rr(ix),drr(ix),rrm(ix)

c----------------------------------------------------------------------|      
c     numerical parameters
      pi = acos(-1.0d0)
      pi4=4.0d0*pi
      pi8=8.0d0*pi
      pi4i=1.0d0/pi4
      pi8i=5.0d-1*pi4i
c
c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do i=1,ix
         rosc(i)=dv(i)*ro(i)
         rxsc(i)=dv(i)*ro(i)*vx(i)
         rysc(i)=dv(i)*ro(i)*vy(i)*rr(i)
         bysc(i)=dv(i)*by(i)/rr(i)
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do i=1,ix-1
         row(i,1)=ro(i)
         vxw(i,1)=vx(i)
         vyw(i,1)=vy(i)
         bxw(i,1)=bxm(i)
         byw(i,1)=by(i)
         row(i,2)=ro(i+1)
         vxw(i,2)=vx(i+1)
         vyw(i,2)=vy(i+1)
         bxw(i,2)=bxm(i)
         byw(i,2)=by(i+1)
      enddo

      call roeflux_mt(fro,frx,fry,fby,cs2,row,vxw,vyw,bxw,byw,ix)

      do i=1,ix-1
        fro(i)=scm(i)*fro(i)
        frx(i)=scm(i)*frx(i)
        fry(i)=scm(i)*fry(i)*rrm(i)
        fby(i)=scm(i)*fby(i)/rrm(i)
      enddo

      do i=2,ix-1
         rosch(i)=rosc(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         rxsch(i)=rxsc(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
         rysch(i)=rysc(i)+0.5d0*dt*( (fry(i-1)-fry(i))/dx(i) )
         bysch(i)=bysc(i)+0.5d0*dt*( (fby(i-1)-fby(i))/dx(i) )
      enddo

      do i=2,ix-1
         b2=bx(i)**2+by(i)**2
         srx=dv(i)
     &   *(ro(i)*gx(i)+(ro(i)*vy(i)**2-by(i)**2*pi4i)/rr(i)*drr(i))
     &       +(cs2*ro(i)+b2*pi8i)*dsc(i)
         rxsch(i)=rxsch(i)+0.5d0*dt*srx
      enddo


c     computation of basic variables on half step

      do i=2,ix-1
         roh(i)=rosch(i)/dv(i)
         vxh(i)=rxsch(i)/rosch(i)
         vyh(i)=rysch(i)/rosch(i)/rr(i)
         byh(i)=bysch(i)/dv(i)*rr(i)
      enddo

c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(vxh,vxw,ix)
      call tvdminmod(vyh,vyw,ix)
      call tvdminmod(byh,byw,ix)

      call roeflux_mt(fro,frx,fry,fby,cs2,row,vxw,vyw,bxw,byw,ix)

      do i=2,ix-2
        fro(i)=scm(i)*fro(i)
        frx(i)=scm(i)*frx(i)
        fry(i)=scm(i)*fry(i)*rrm(i)
        fby(i)=scm(i)*fby(i)/rrm(i)
      enddo

      do i=2,ix-2
         rosc(i)=rosc(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         rxsc(i)=rxsc(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
         rysc(i)=rysc(i)+dt*( (fry(i-1)-fry(i))/dx(i) )
         bysc(i)=bysc(i)+dt*( (fby(i-1)-fby(i))/dx(i) )
      enddo

      do i=3,ix-2
         b2=bx(i)**2+byh(i)**2
         srx=dv(i)
     &   *(roh(i)*gx(i)+(roh(i)*vyh(i)**2-byh(i)**2*pi4i)/rrm(i)*drr(i))
     &       +(cs2*roh(i)+b2*pi8i)*dsc(i)
         rxsc(i)=rxsc(i)+dt*srx
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         ro(i)=rosc(i)/dv(i)
         vx(i)=rxsc(i)/rosc(i)
         vy(i)=rysc(i)/rosc(i)/rr(i)
         by(i)=bysc(i)/dv(i)*rr(i)
      enddo
c----------------------------------------------------------------------|

      return
      end
