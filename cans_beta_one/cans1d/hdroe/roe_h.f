c======================================================================|
      subroutine roe_h(ro,pr,vx,dt,gm,dx,ix)
c======================================================================|
c
c NAME  roe_h
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * hydrodynamics
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
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

      dimension ro(ix),pr(ix),vx(ix)
      dimension ee(ix),rx(ix)

      dimension fro(ix),fee(ix),frx(ix)
      dimension roh(ix),eeh(ix),rxh(ix)
      dimension prh(ix),vxh(ix)
      dimension row(ix,2),prw(ix,2),vxw(ix,2)

c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do i=1,ix
         rx(i)=ro(i)*vx(i)
         v2=vx(i)**2
         ee(i)=pr(i)/(gm-1.0d0) +0.5d0*ro(i)*v2
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do i=1,ix-1
         row(i,1)=ro(i)
         prw(i,1)=pr(i)
         vxw(i,1)=vx(i)
         row(i,2)=ro(i+1)
         prw(i,2)=pr(i+1)
         vxw(i,2)=vx(i+1)
      enddo

      call roeflux_h(fro,fee,frx,gm,row,prw,vxw,ix)

      do i=2,ix-1
         roh(i)=ro(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         eeh(i)=ee(i)+0.5d0*dt*( (fee(i-1)-fee(i))/dx(i) )
         rxh(i)=rx(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo

c     computation of basic variables on half step

      do i=2,ix-1
         vxh(i)=rxh(i)/roh(i)
         v2=vxh(i)**2
         prh(i)=(gm-1.0d0)*(eeh(i)-0.5d0*roh(i)*v2)
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(prh,prw,ix)
      call tvdminmod(vxh,vxw,ix)

      call roeflux_h(fro,fee,frx,gm,row,prw,vxw,ix)

      do i=3,ix-2
         ro(i)=ro(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         ee(i)=ee(i)+dt*( (fee(i-1)-fee(i))/dx(i) )
         rx(i)=rx(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo
c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         vx(i)=rx(i)/ro(i)
         v2=vx(i)**2
         pr(i)=(gm-1.0d0)*(ee(i)-0.5d0*ro(i)*v2)
      enddo
c----------------------------------------------------------------------|
      return
      end
