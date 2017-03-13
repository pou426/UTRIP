c======================================================================|
      subroutine mlw_h_o(ro,pr,vx,dt,qav,gm
     &     ,hh,hhm,gi,gim,gd,gdm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_h_o
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * hydrodynamics
c        * orthogonal general coordinate
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
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
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c    hh(ix,3), hhm(ix,3) : [double] metric, h1,h2,h3
c    gi(ix,3), gim(ix,3) : [double] curvature, G23, G31, G12
c    gd(ix,3), gdm(ix,3) : [double] curvature, G32, G13, G21
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension roh(ix),prh(ix),vxh(ix)
      dimension rosc(ix),eesc(ix),rxsc(ix)
      dimension rosch(ix),eesch(ix),rxsch(ix)
      dimension drosc(ix),deesc(ix),drxsc(ix)
      dimension fx(ix),qx(ix)
      dimension ss(ix)
      dimension hh(ix,3),gi(ix,3),gd(ix,3)
      dimension hhm(ix,3),gim(ix,3),gdm(ix,3)

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
         rosc(i) = 0.d0
         eesc(i) = 0.d0
         rxsc(i) = 0.d0
         drosc(i) = 0.d0
         deesc(i) = 0.d0
         drxsc(i) = 0.d0
         rosch(i) = 0.d0
         eesch(i) = 0.d0
         rxsch(i) = 0.d0
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         vv=vx(i)**2
         ee    = pr(i)/(gm-1)+0.5*ro(i)*vv
         rx    = ro(i)*vx(i)

         sc=hh(i,1)*hh(i,2)*hh(i,3)
         rosc(i) = ro(i)*sc
         eesc(i) = ee*sc
         rxsc(i) = rx*sc
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
      enddo
      call mlwhalf(rosc ,rosch ,drosc,dt,fx,dxi,dxim,ix)

c---  energy ---
      do i=1,ix
         vv=vx(i)**2
         ep=pr(i)*gm/(gm-1.)+0.5*ro(i)*vv
         fx(i)= ep*vx(i)
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
      enddo
      call mlwhalf(eesc ,eesch ,deesc ,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix
         t11= ro(i)*vx(i)**2+pr(i)
         t22=pr(i)
         t33=pr(i)
         t12=0.d0
         t13=0.d0
         fx(i)= t11
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
         ss(i)=-gd(i,3)*t22-gi(i,2)*t33+gi(i,3)*t12+gd(i,2)*t13
         ss(i)= ss(i)*hh(i,1)*hh(i,2)*hh(i,3)
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix)
      call mlwsrch(rxsch,drxsc,dt,ss,ix)
c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do i=1,ix-1
         scm=hhm(i,1)*hhm(i,2)*hhm(i,3)
         roh(i) = rosch(i)/scm
         eeh    = eesch(i)/scm
         rxh    = rxsch(i)/scm
         vxh(i) = rxh/roh(i)
         vv=vxh(i)**2
         prh(i)   = (gm-1)*(eeh-0.5d0*roh(i)*vv)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)
         fx(i)= fx(i)*hhm(i,2)*hhm(i,3)
      enddo
      call mlwfull(drosc ,dt,fx,dxi,ix)

c---  energy     ---
      do i=1,ix-1
         vv=vxh(i)**2
         eph   = prh(i)*gm/(gm-1.)+0.5*roh(i)*vv
         fx(i)= eph*vxh(i)
         fx(i)= fx(i)*hhm(i,2)*hhm(i,3)
      enddo
      call mlwfull(deesc ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         t11=prh(i)+ roh(i)*vxh(i)**2
         t22=prh(i)
         t33=prh(i)
         t12=0.d0
         t13=0.d0
         fx(i)= t11
         fx(i)= fx(i)*hhm(i,2)*hhm(i,3)
         ss(i)=-gdm(i,3)*t22-gim(i,2)*t33+gim(i,3)*t12+gdm(i,2)*t13
         ss(i)= ss(i)*hhm(i,1)*hhm(i,2)*hhm(i,3)
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
         rosc(i) = rosc(i) +drosc(i)
         eesc(i) = eesc(i) +deesc(i)
         rxsc(i) = rxsc(i) +drxsc(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         sc=hh(i,1)*hh(i,2)*hh(i,3)
         ro(i)=rosc(i)/sc
         ee   =eesc(i)/sc
         rx   =rxsc(i)/sc
         vx(i) = rx/ro(i)
         vv=vx(i)**2
         pr(i) = (gm-1)*(ee - 0.5*ro(i)*vv)
      enddo

      return
      end
