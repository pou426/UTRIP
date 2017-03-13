c======================================================================|
      subroutine roe_m3t(ro,vx,vy,vz,by,bz,bx,bxm,dt,cs2,dx,ix)
c======================================================================|
c
c NAME  roe_m3t
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * isothermal 3-component MHD
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    vz(ix): [double] velocity 
c    by(ix): [double] magnetic field
c    bz(ix): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    cs2: [double] square of sound speed
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)
      dimension ro(ix),vx(ix),vy(ix),vz(ix),by(ix),bz(ix)
      dimension rx(ix),ry(ix),rz(ix)

      dimension bx(ix),bxm(ix)

      dimension fro(ix),frx(ix),fry(ix),frz(ix),fby(ix),fbz(ix)
      dimension roh(ix),rxh(ix),ryh(ix),rzh(ix),byh(ix),bzh(ix)
      dimension vxh(ix),vyh(ix),vzh(ix)

      dimension row(ix,2),vxw(ix,2),vyw(ix,2),vzw(ix,2)
     &         ,bxw(ix,2),byw(ix,2),bzw(ix,2)

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
         rz(i)=ro(i)*vz(i)
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
         vzw(i,1)=vz(i)
         bzw(i,1)=bz(i)
         row(i,2)=ro(i+1)
         vxw(i,2)=vx(i+1)
         vyw(i,2)=vy(i+1)
         bxw(i,2)=bxm(i)
         byw(i,2)=by(i+1)
         vzw(i,2)=vz(i+1)
         bzw(i,2)=bz(i+1)
      enddo

      call roeflux_m3t(fro,frx,fry,frz,fby,fbz,cs2
     &              ,row,vxw,vyw,vzw,bxw,byw,bzw,ix)

      do i=2,ix-1
         roh(i)=ro(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         rxh(i)=rx(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
         ryh(i)=ry(i)+0.5d0*dt*( (fry(i-1)-fry(i))/dx(i) )
         byh(i)=by(i)+0.5d0*dt*( (fby(i-1)-fby(i))/dx(i) )
         rzh(i)=rz(i)+0.5d0*dt*( (frz(i-1)-frz(i))/dx(i) )
         bzh(i)=bz(i)+0.5d0*dt*( (fbz(i-1)-fbz(i))/dx(i) )
      enddo

c     computation of basic variables on half step

      do i=2,ix-1
         vxh(i)=rxh(i)/roh(i)
         vyh(i)=ryh(i)/roh(i)
         vzh(i)=rzh(i)/roh(i)
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(vxh,vxw,ix)
      call tvdminmod(vyh,vyw,ix)
      call tvdminmod(byh,byw,ix)
      call tvdminmod(vzh,vzw,ix)
      call tvdminmod(bzh,bzw,ix)

      call roeflux_m3t(fro,frx,fry,frz,fby,fbz,cs2
     &              ,row,vxw,vyw,vzw,bxw,byw,bzw,ix)

      do i=3,ix-2
         ro(i)=ro(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         rx(i)=rx(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
         ry(i)=ry(i)+dt*( (fry(i-1)-fry(i))/dx(i) )
         by(i)=by(i)+dt*( (fby(i-1)-fby(i))/dx(i) )
         rz(i)=rz(i)+dt*( (frz(i-1)-frz(i))/dx(i) )
         bz(i)=bz(i)+dt*( (fbz(i-1)-fbz(i))/dx(i) )
      enddo
c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         vx(i)=rx(i)/ro(i)
         vy(i)=ry(i)/ro(i)
         vz(i)=rz(i)/ro(i)
      enddo
c----------------------------------------------------------------------|
      return
      end
