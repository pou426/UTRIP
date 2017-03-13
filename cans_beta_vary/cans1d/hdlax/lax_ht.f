c======================================================================|
      subroutine lax_ht(ro,vx,dt,cs2,dx,dxm,ix)
c======================================================================|
c
c NAME  lax_ht
c
c PURPOSE
c    solve eqs. by Lax-Friedrich method with effects of
c        * isothermal hydrodynamics
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity along the x-cordinate
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    cs2: [double] square of sound speed
c    ix: [integer] dimension size
c
c HISTORY
c    written 2006-8-2 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension rx(ix)
      dimension ron(ix),rxn(ix)
      dimension fx(ix)

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         pr(i) = cs2*ro(i)
         rx(i) = ro(i)*vx(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)
      enddo
      call lax(ron,ro,dt,fx,dx,dxm,ix)

c---  x-momentum ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)**2+pr(i)
      enddo
      call lax(rxn,rx,dt,fx,dx,dxm,ix)

c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = ron(i)
         rx(i) = rxn(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
      enddo

      return
      end
