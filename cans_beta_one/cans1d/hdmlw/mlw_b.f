c======================================================================|
      subroutine mlw_b(u,dt,qav,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_b
c
c PURPOSE
c    solve Burgers eq. by modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    u(ix): [double] 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    dxm(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2006-8-2 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension u(ix)
      dimension uh(ix)
      dimension du(ix)
      dimension fx(ix),qx(ix)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.0/dxm(i)
         dxi(i) = 1.0/dx(i)
      enddo

c----------------------------------------------------------------------|
c     initialize du etc.                                   
c----------------------------------------------------------------------|
      do i=1,ix
         du(i) = 0.0
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
      do i=1,ix
         fx(i)= u(i)*u(i)/2
      enddo
      call mlwhalf(u ,uh ,du,dt,fx,dxi,dxim,ix)

c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= uh(i)*uh(i)/2
      enddo
      call mlwfull(du ,dt,fx,dxi,ix)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0e0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(u(i+1)-u(i))-1.0e-4)
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(ro,du,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         u(i)= u(i) +du(i)
      enddo

      return
      end
