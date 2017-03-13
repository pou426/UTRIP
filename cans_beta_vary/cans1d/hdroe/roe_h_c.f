c======================================================================|
      subroutine roe_h_c(ro,pr,vx,dt,gm,dsc,scm,dv,dx,ix)
c======================================================================|
c
c NAME  roe_h_c
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * hydrodynamics
c        * non-uniform cross section
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
c    scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
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

      dimension dsc(ix),scm(ix),dv(ix)

      dimension rosc(ix),eesc(ix),rxsc(ix)
      dimension rosch(ix),eesch(ix),rxsch(ix)

      dimension fro(ix),fee(ix),frx(ix)
      dimension roh(ix)
      dimension prh(ix),vxh(ix)
      dimension row(ix,2),prw(ix,2),vxw(ix,2)


c----------------------------------------------------------------------|
c     computation of conservative variables

      do i=1,ix
         rosc(i)=dv(i)*ro(i)
         rxsc(i)=dv(i)*ro(i)*vx(i)
         v2=vx(i)**2
         eesc(i)=dv(i)*(pr(i)/(gm-1.0d0)+0.5d0*ro(i)*v2)
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
      do i=1,ix-1
        fro(i)=scm(i)*fro(i)
        fee(i)=scm(i)*fee(i)
        frx(i)=scm(i)*frx(i)
      enddo

      do i=2,ix-1
         rosch(i)=rosc(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         eesch(i)=eesc(i)+0.5d0*dt*( (fee(i-1)-fee(i))/dx(i) )
         rxsch(i)=rxsc(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo

      do i=2,ix-1
         srx=pr(i)*dsc(i)
         rxsch(i)=rxsch(i)+0.5d0*dt*srx
      enddo

c     computation of basic variables on half step

      do i=2,ix-1
         roh(i)=rosch(i)/dv(i)
         vxh(i)=rxsch(i)/rosch(i)
         v2=vxh(i)**2
         prh(i)=(gm-1.0d0)* (eesch(i)/dv(i)-0.5d0*roh(i)*v2)
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(prh,prw,ix)
      call tvdminmod(vxh,vxw,ix)

      call roeflux_h(fro,fee,frx,gm,row,prw,vxw,ix)
      do i=2,ix-2
        fro(i)=scm(i)*fro(i)
        fee(i)=scm(i)*fee(i)
        frx(i)=scm(i)*frx(i)
      enddo

      do i=3,ix-2
         rosc(i)=rosc(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         eesc(i)=eesc(i)+dt*( (fee(i-1)-fee(i))/dx(i) )
         rxsc(i)=rxsc(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo

      do i=3,ix-2
        srx=prh(i)*dsc(i)
        rxsc(i)=rxsc(i)+dt*srx
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         ro(i)=rosc(i)/dv(i)
         vx(i)=rxsc(i)/rosc(i)
         v2=vx(i)**2
         pr(i)=(gm-1.0d0)* (eesc(i)/dv(i)-0.5d0*ro(i)*v2)
      enddo
c----------------------------------------------------------------------|
      return
      end
