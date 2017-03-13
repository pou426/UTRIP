c======================================================================|
      subroutine roe_m_g(ro,pr,vx,vy,by,bx,bxm,dt,gm,gx,gy,dx,ix)
c======================================================================|
c
c NAME  roe_m_g
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * MHD
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
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
c    gx(ix) : [double] gravity
c    gy(ix) : [double] gravity
c    dx(ix) : [double] grid spacing
c    gm: [double] polytropic index gamma
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)

      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
      dimension ee(ix),rx(ix),ry(ix)

      dimension bx(ix),bxm(ix)

      dimension fro(ix),fee(ix),frx(ix),fry(ix),fby(ix)
      dimension roh(ix),eeh(ix),rxh(ix),ryh(ix),byh(ix)
      dimension prh(ix),vxh(ix),vyh(ix)
      dimension row(ix,2),prw(ix,2),vxw(ix,2),vyw(ix,2)
     &         ,bxw(ix,2),byw(ix,2)

      dimension gx(ix),gy(ix)

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
         rx(i)=ro(i)*vx(i)
         ry(i)=ro(i)*vy(i)
         v2=vx(i)**2+vy(i)**2
         b2=bx(i)**2+by(i)**2
         ee(i)=pr(i)/(gm-1.0d0) +0.5d0*ro(i)*v2 + pi8i*b2
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do i=1,ix-1
         row(i,1)=ro(i)
         prw(i,1)=pr(i)
         vxw(i,1)=vx(i)
         vyw(i,1)=vy(i)
         bxw(i,1)=bxm(i)
         byw(i,1)=by(i)
         row(i,2)=ro(i+1)
         prw(i,2)=pr(i+1)
         vxw(i,2)=vx(i+1)
         vyw(i,2)=vy(i+1)
         bxw(i,2)=bxm(i)
         byw(i,2)=by(i+1)
      enddo

      call roeflux_m(fro,fee,frx,fry,fby,gm,row,prw,vxw,vyw,bxw,byw,ix)

      do i=2,ix-1
         roh(i)=ro(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         eeh(i)=ee(i)+0.5d0*dt*( (fee(i-1)-fee(i))/dx(i) )
         rxh(i)=rx(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
         ryh(i)=ry(i)+0.5d0*dt*( (fry(i-1)-fry(i))/dx(i) )
         byh(i)=by(i)+0.5d0*dt*( (fby(i-1)-fby(i))/dx(i) )
      enddo

      do i=2,ix-1
         see=ro(i)*(vx(i)*gx(i)+vy(i)*gy(i))
         eeh(i)=eeh(i)+0.5d0*dt*see
         srx=ro(i)*gx(i)
         rxh(i)=rxh(i)+0.5d0*dt*srx
         sry=ro(i)*gy(i)
         ryh(i)=ryh(i)+0.5d0*dt*sry
      enddo

c     computation of basic variables on half step

      do i=2,ix-1
         vxh(i)=rxh(i)/roh(i)
         vyh(i)=ryh(i)/roh(i)
         v2=vxh(i)**2+vyh(i)**2
         b2= bx(i)**2+byh(i)**2
         prh(i)=(gm-1.0d0)*(eeh(i)-0.5d0*roh(i)*v2 -pi8i*b2)
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(prh,prw,ix)
      call tvdminmod(vxh,vxw,ix)
      call tvdminmod(vyh,vyw,ix)
      call tvdminmod(byh,byw,ix)

      call roeflux_m(fro,fee,frx,fry,fby,gm,row,prw,vxw,vyw,bxw,byw,ix)

      do i=3,ix-2
         ro(i)=ro(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         ee(i)=ee(i)+dt*( (fee(i-1)-fee(i))/dx(i) )
         rx(i)=rx(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
         ry(i)=ry(i)+dt*( (fry(i-1)-fry(i))/dx(i) )
         by(i)=by(i)+dt*( (fby(i-1)-fby(i))/dx(i) )
      enddo

      do i=3,ix-2
         see=roh(i)*(vxh(i)*gx(i)+vyh(i)*gy(i))
         ee(i)=ee(i)+dt*see
         srx=roh(i)*gx(i)
         rx(i)=rx(i)+dt*srx
         sry=roh(i)*gy(i)
         ry(i)=ry(i)+dt*sry
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         vx(i)=rx(i)/ro(i)
         vy(i)=ry(i)/ro(i)
         v2=vx(i)**2+vy(i)**2
         b2=bx(i)**2+by(i)**2
         pr(i)=(gm-1.0d0)*(ee(i)-0.5d0*ro(i)*v2 -pi8i*b2)
      enddo
c----------------------------------------------------------------------|
      return
      end
