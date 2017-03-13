c======================================================================|
      subroutine mlw_h_cg(ro,pr,vx,dt,qav,gm,gx,gxm,
     &             sc,dsc,scm,dscm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_h_cg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * hydrodynamics
c        * non-uniform cross section
c        * gravity
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
c    gx(ix), gxm(ix) : [double] gravity
c    dx(ix), dxm(ix): [double] grid spacing
c    sc(ix), scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
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

      dimension ro(ix),pr(ix),vx(ix)
      dimension ee(ix),rx(ix)

      dimension roh(ix),eeh(ix)
      dimension prh(ix),vxh(ix)

      dimension sc(ix),scm(ix)
      dimension dsc(ix),dscm(ix)
      dimension rosc(ix),eesc(ix),rxsc(ix)
      dimension rosch(ix),eesch(ix),rxsch(ix)
      dimension drosc(ix),deesc(ix),drxsc(ix)

      dimension fx(ix),qx(ix)
      dimension ss(ix)

      dimension gx(ix), gxm(ix)

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
         drxsc(i )= 0.0
         deesc(i) = 0.0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix       
         vv=vx(i)**2
         ee(i) = pr(i)/(gm-1)+0.5*ro(i)*vv
         rx(i) = ro(i)*vx(i)
      enddo
      do i=1,ix       
         rosc(i) = ro(i)*sc(i)
         rxsc(i) = rx(i)*sc(i)
         eesc(i) = ee(i)*sc(i)
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

c---  energy ---
      do i=1,ix       
c        vv=vx(i)**2
c        ep=pr(i)*gm/(gm-1.)+0.5*ro(i)*vv
         ep=pr(i)+ee(i)
         fx(i)= ep*vx(i)
         ss(i)= ro(i)*vx(i)*gx(i)
         fx(i)= fx(i)*sc(i)
         ss(i)= ss(i)*sc(i)
      enddo
      call mlwhalf(eesc ,eesch ,deesc ,dt,fx,dxi,dxim,ix)
      call mlwsrch(eesch ,deesc ,dt,ss,ix)

c---  x-momentum ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)**2+pr(i)
         ss(i)= ro(i)*gx(i)
         fx(i)= fx(i)*sc(i)
         ss(i)= ss(i)*sc(i)+pr(i)*dsc(i)
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix)
      call mlwsrch(rxsch,drxsc,dt,ss,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         roh(i)   = rosch(i)/scm(i)
         eeh(i)   = eesch(i)/scm(i)
         vxh(i)   = rxsch(i)/rosch(i)
         vv=vxh(i)**2
         prh(i)   = (gm-1)*(eeh(i)-0.5*roh(i)*vv)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)
         fx(i)= fx(i)*scm(i)
      enddo
      call mlwfull(drosc ,dt,fx,dxi,ix)

c---  energy     ---
      do i=1,ix-1
c        vv=vxh(i)**2
c        eph   = prh(i)*gm/(gm-1.)+0.5*roh(i)*vv
         eph   = prh(i)+eeh(i)
         fx(i)= eph*vxh(i)
         fx(i)= fx(i)*scm(i)
         ss(i)= roh(i)*vxh(i)*gxm(i)
         ss(i)= ss(i)*scm(i)
      enddo
      call mlwfull(deesc ,dt,fx,dxi,ix)
      call mlwsrcf(deesc ,dt,ss,ux0,ux1,ix)

c---  x-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)**2+prh(i)
         ss(i)= roh(i)*gxm(i)
         fx(i)= fx(i)*scm(i)
         ss(i)= ss(i)*scm(i)+prh(i)*dscm(i)
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ix)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix)
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
      call mlwartv(rosc,drosc,dt,qx,dxi,dxim,ix)
      call mlwartv(eesc,deesc,dt,qx,dxi,dxim,ix)
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = (rosc(i) +drosc(i))/sc(i)
         rx(i) = (rxsc(i) +drxsc(i))/sc(i)
         ee(i) = (eesc(i) +deesc(i))/sc(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
         vv=vx(i)**2
         pr(i) = (gm-1)*(ee(i) - 0.5*ro(i)*vv)
      enddo

      return
      end
