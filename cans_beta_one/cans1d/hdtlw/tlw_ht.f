c======================================================================|
      subroutine tlw_ht(ro,vx,dt,qav,cs2,dx,dxm,ix)
c======================================================================|
c
c NAME  tlw_h
c
c PURPOSE
c    solve eqs. by two-step Lax-Wendroff method with effects of
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
c    written 2002-8-2 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension rx(ix)
      dimension roh(ix),rxh(ix)
      dimension prh(ix),vxh(ix)
      dimension dro(ix),drx(ix)
      dimension fx(ix),qx(ix)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.0/dxm(i)
         dxi(i) = 1.0/dx(i)
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do i=1,ix
         dro(i) = 0.0
         drx(i) = 0.0
      enddo

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
      call tlwhalf(ro ,roh ,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)**2+pr(i)
      enddo
      call tlwhalf(rx,rxh,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do i=1,ix-1
         vxh(i)   = rxh(i)/roh(i)
         prh(i)   = cs2*roh(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)
      enddo
      call tlwfull(dro ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)**2+prh(i)
      enddo
      call tlwfull(drx,dt,fx,dxi,ix)
c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
      zero=0.0e0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call tlwartv(ro,dro,dt,qx,dxi,dxim,ix)
      call tlwartv(rx,drx,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = ro(i) +dro(i)
         rx(i) = rx(i) +drx(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
      enddo

      return
      end
