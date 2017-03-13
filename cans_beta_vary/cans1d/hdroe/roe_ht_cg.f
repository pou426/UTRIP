c======================================================================|
      subroutine roe_ht_cg(ro,vx,dt,cs2,gx,dsc,scm,dv,dx,ix)
c======================================================================|
c
c NAME  roe_ht_cg
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * isothermal hydrodynamics
c        * non-uniform cross section
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    gx(ix), gxm(ix) : [double] gravity
c    scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    dx(ix) : [double] grid spacing
c    cs2: [double] square of sound speed
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on N. Fukuda's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix)

      dimension ro(ix),vx(ix)

      dimension dsc(ix),scm(ix),dv(ix)

      dimension rosc(ix),rxsc(ix)
      dimension rosch(ix),rxsch(ix)

      dimension fro(ix),frx(ix)
      dimension roh(ix),vxh(ix)
      dimension row(ix,2),vxw(ix,2)

      dimension gx(ix)

c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do i=1,ix
         rosc(i)=dv(i)*ro(i)
         rxsc(i)=dv(i)*ro(i)*vx(i)
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do i=1,ix-1
         row(i,1)=ro(i)
         vxw(i,1)=vx(i)
         row(i,2)=ro(i+1)
         vxw(i,2)=vx(i+1)
      enddo

      call roeflux_ht(fro,frx,cs2,row,vxw,ix)
      do i=1,ix-1
        fro(i)=scm(i)*fro(i)
        frx(i)=scm(i)*frx(i)
      enddo

      do i=2,ix-1
         rosch(i)=rosc(i)+0.5d0*dt*( (fro(i-1)-fro(i))/dx(i) )
         rxsch(i)=rxsc(i)+0.5d0*dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo

      do i=2,ix-1
         srx=dv(i)*ro(i)*gx(i)+cs2*ro(i)*dsc(i)
         rxsch(i)=rxsch(i)+0.5d0*dt*srx
      enddo

c     computation of basic variables on half step

      do i=2,ix-1
         roh(i)=rosch(i)/dv(i)
         vxh(i)=rxsch(i)/rosch(i)
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
      call tvdminmod(roh,row,ix)
      call tvdminmod(vxh,vxw,ix)

      call roeflux_ht(fro,frx,cs2,row,vxw,ix)
      do i=2,ix-2
        fro(i)=scm(i)*fro(i)
        frx(i)=scm(i)*frx(i)
      enddo

      do i=3,ix-2
         rosc(i)=rosc(i)+dt*( (fro(i-1)-fro(i))/dx(i) )
         rxsc(i)=rxsc(i)+dt*( (frx(i-1)-frx(i))/dx(i) )
      enddo

      do i=3,ix-2
         srx=dv(i)*roh(i)*gx(i)+cs2*roh(i)*dsc(i)
         rxsc(i)=rxsc(i)+dt*srx
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do i=3,ix-2
         ro(i)=rosc(i)/dv(i)
         vx(i)=rxsc(i)/rosc(i)
      enddo
c----------------------------------------------------------------------|
      return
      end
