c======================================================================|
      subroutine mlw_ht_c(ro,vx,dt,qav,cs2,sc,dsc,scm,dscm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_ht_c
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal hydrodynamics
c        * non-uniform cross section
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
c    sc(ix), scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    cs2: [double] square of sound speed
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ux0(ix),ux1(ix)

      dimension ro(ix),vx(ix)
      dimension pr(ix),rx(ix)

      dimension roh(ix)
      dimension prh(ix),vxh(ix)

      dimension sc(ix),scm(ix)
      dimension dsc(ix),dscm(ix)
      dimension rosc(ix),rxsc(ix)
      dimension rosch(ix),rxsch(ix)
      dimension drosc(ix),drxsc(ix)

      dimension fx(ix),qx(ix)
      dimension ss(ix)
c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.0/dxm(i)
         dxi(i) = 1.0/dx(i)
      enddo
      do i=2,ix-1
         ux1(i)  = 0.5*dxm(i-1)/dx(i)
         ux0(i)  = 0.5*dxm(i)/dx(i)
      enddo
c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      
      do i=1,ix
         drosc(i) = 0.0
         drxsc(i)= 0.0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix       
         rx(i) = ro(i)*vx(i)
         pr(i) = cs2*ro(i)
      enddo
c----------------------------------------------------------------------|
c     
c----------------------------------------------------------------------|
      do i=1,ix       
         rosc(i) = ro(i)*sc(i)
         rxsc(i) = rx(i)*sc(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)
         fx(i)= fx(i)*sc(i)
      enddo
      call mlwhalf(rosc ,rosch ,drosc ,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)**2+pr(i)
         fx(i)= fx(i)*sc(i)
         ss(i)= pr(i)*dsc(i)
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix)
      call mlwsrch(rxsch,drxsc,dt,ss,ix)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         roh(i)   = rosch(i)/scm(i)
         vxh(i)   = rxsch(i)/rosch(i)
         prh(i)   = cs2*roh(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)
         fx(i)= scm(i)*fx(i)
      enddo
      call mlwfull(drosc ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)**2+prh(i)
         fx(i)= scm(i)*fx(i)
         ss(i)= prh(i)*dscm(i)
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ix)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix)

c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0e0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
      call mlwartv(rosc,drosc,dt,qx,dxi,dxim,ix)
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i)= (rosc(i) +drosc(i))/sc(i)
         rx(i)= (rxsc(i) +drxsc(i))/sc(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
      enddo

      return
      end
