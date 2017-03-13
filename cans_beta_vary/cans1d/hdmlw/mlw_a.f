c======================================================================|
      subroutine mlw_a(ro,dt,qav,vx,vxm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_a
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * simple advection
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vx(ix), vxm(ix) : [double] velocity along the x-cordinate
c    dxm(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ro(ix),vx(ix)
      dimension roh(ix)
      dimension vxm(ix)
      dimension dro(ix)
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
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)
      enddo
      call mlwhalf(ro ,roh ,dro,dt,fx,dxi,dxim,ix)

c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxm(i)
      enddo
      call mlwfull(dro ,dt,fx,dxi,ix)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0e0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(ro,dro,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i)= ro(i) +dro(i)
      enddo

      return
      end
