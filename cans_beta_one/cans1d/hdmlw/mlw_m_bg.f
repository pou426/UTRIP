c======================================================================|
      subroutine mlw_m_bg(ro,pr,vx,vy,by,bx,bxm,dt,qav,gm
     &             ,gx,gxm,sc,dsc,scm,dscm,rr,rrm,drr,drrm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_m_bg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * axial symmetry
c        * non-uniform poloidal magnetic field
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity along the x-cordinate
c    vy(ix): [double] velocity in the rotating direction
c    by(ix): [double] magnetic field in rotation direction
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    gx(ix), gxm(ix) : [double] gravity
c    sc(ix), scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    rr(ix), rrm(ix) : [double] distance from rotation axis
c    drr(ix), drrm(ix) : [double] distance gradient from rotation axis
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dxm(ix),dx(ix)
      dimension dxim(ix),dxi(ix)
      dimension ux0(ix),ux1(ix)

      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
      dimension ee(ix),rx(ix),ry(ix),ez(ix)

      dimension roh(ix),prh(ix),vxh(ix),vyh(ix),byh(ix)
      dimension eeh(ix),ezh(ix)

      dimension bx(ix),bxm(ix)

      dimension sc(ix),scm(ix),dsc(ix),dscm(ix)
      dimension rosc(ix),eesc(ix),rxsc(ix),rysc(ix),bysc(ix)
      dimension rosch(ix),eesch(ix),rxsch(ix),rysch(ix),bysch(ix)
      dimension drosc(ix),deesc(ix),drxsc(ix),drysc(ix),dbysc(ix)

      dimension fx(ix),qx(ix)
      dimension ss(ix)

      dimension gx(ix),gxm(ix)
      dimension rr(ix),rrm(ix),drr(ix),drrm(ix)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
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
         deesc(i) = 0.0
         drxsc(i) = 0.0
         drysc(i) = 0.0
         dbysc(i) = 0.0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         vv=vx(i)**2+vy(i)**2
c         bb=by(i)**2
         bb=by(i)**2+bx(i)**2
         ee(i) = pr(i)/(gm-1)+0.5*ro(i)*vv+pi8i*bb
         rx(i) = ro(i)*vx(i)
         ry(i) = ro(i)*vy(i)
      enddo
      do i=1,ix
         ez(i) = -vx(i)*by(i)+vy(i)*bx(i)
      enddo
      do i=1,ix       
         rosc(i) = ro(i)*sc(i)
         eesc(i) = ee(i)*sc(i)
         rxsc(i) = rx(i)*sc(i)
         rysc(i) = ry(i)*sc(i)*rr(i)
         bysc(i) = by(i)*sc(i)/rr(i)
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
         vv=vx(i)**2+vy(i)**2
         ep    = pr(i)*gm/(gm-1.)+0.5*ro(i)*vv
         fx(i)= ep*vx(i)-by(i)*ez(i)*pi4i
         fx(i)= fx(i)*sc(i)
         ss(i)= ro(i)*vx(i)*gx(i)
         ss(i)= ss(i)*sc(i)
      enddo
      call mlwhalf(eesc ,eesch ,deesc ,dt,fx,dxi,dxim,ix)
      call mlwsrch(eesch ,deesc ,dt,ss,ix)

c---  x-momentum ---
      do i=1,ix
c         bbm=by(i)**2
         bbm=by(i)**2-bx(i)**2
c         bb =by(i)**2
         bb =by(i)**2+bx(i)**2
         fx(i)= ro(i)*vx(i)**2+pr(i)+bbm*pi8i
         ss(i)= ro(i)*gx(i)
     &         +(ro(i)*vy(i)**2-by(i)**2*pi4i)/rr(i)*drr(i)
         fx(i)= fx(i)*sc(i)
         ss(i)= ss(i)*sc(i)+(pr(i)+bb*pi8i)*dsc(i)
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix)
      call mlwsrch(rxsch,drxsc,dt,ss,ix)

c---  y-momentum ---
      do i=1,ix
         fx(i)= ro(i)*vx(i)*vy(i)-bx(i)*by(i)*pi4i
         fx(i)= fx(i)*sc(i)*rr(i)
      enddo
      call mlwhalf(rysc,rysch,drysc,dt,fx,dxi,dxim,ix)

c---  y-magnetic ---
      do i=1,ix
         fx(i)= -ez(i)
         fx(i)= fx(i)*sc(i)
         fx(i)= fx(i)/rr(i)
      enddo
      call mlwhalf(bysc ,bysch ,dbysc ,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         roh(i)   = rosch(i)/scm(i)
         eeh(i)   = eesch(i)/scm(i)
         vxh(i)   = rxsch(i)/rosch(i)
         vyh(i)   = rysch(i)/rosch(i)/rrm(i)
         byh(i)   = bysch(i)/scm(i)*rrm(i)
         vv=vxh(i)**2+vyh(i)**2
c         bb=byh(i)**2
         bb=bxm(i)**2+byh(i)**2
         prh(i) = (gm-1)*(eeh(i)-0.5*roh(i)*vv-pi8i*bb)
      enddo
      do i=1,ix
        ezh(i)=-vxh(i)*byh(i)+vyh(i)*bxm(i)
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
         vv=vxh(i)**2+vyh(i)**2
         eph = prh(i)*gm/(gm-1.)+0.5*roh(i)*vv
         fx(i)= eph*vxh(i)-byh(i)*ezh(i)*pi4i
         ss(i)= roh(i)*vxh(i)*gxm(i)
         fx(i)= fx(i)*scm(i)
         ss(i)= ss(i)*scm(i)
      enddo
      call mlwfull(deesc ,dt,fx,dxi,ix)
      call mlwsrcf(deesc ,dt,ss,ux0,ux1,ix)

c---  x-momentum ---
      do i=1,ix-1
c         bbm=byh(i)**2
c         bb =byh(i)**2
         bbm=byh(i)**2-bxm(i)**2
         bb =byh(i)**2+bxm(i)**2
         fx(i)= roh(i)*vxh(i)**2+prh(i)+bbm*pi8i
         ss(i)= roh(i)*gxm(i)
     &        +(roh(i)*vyh(i)**2-byh(i)**2*pi4i)/rrm(i)*drrm(i)
         fx(i)= fx(i)*scm(i)
         ss(i)= ss(i)*scm(i)+(prh(i)+bb*pi8i)*dscm(i)
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ix)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix)

c---  y-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)*vyh(i)-bxm(i)*byh(i)*pi4i
         fx(i)= fx(i)*scm(i)*rrm(i)
      enddo
      call mlwfull(drysc,dt,fx,dxi,ix)

c---  y-magnetic ---
      do i=1,ix-1
         fx(i)= -ezh(i)
         fx(i)= fx(i)*scm(i)
         fx(i)= fx(i)/rrm(i)
      enddo
      call mlwfull(dbysc ,dt,fx,dxi,ix)
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0e0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
      call mlwartv(rosc,drosc,dt,qx,dxi,dxim,ix)
      call mlwartv(eesc,deesc,dt,qx,dxi,dxim,ix)
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix)
      call mlwartv(rysc,drysc,dt,qx,dxi,dxim,ix)
      call mlwartv(bysc,dbysc,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = (rosc(i) +drosc(i))/sc(i)
         ee(i) = (eesc(i) +deesc(i))/sc(i)
         rx(i) = (rxsc(i) +drxsc(i))/sc(i)
         ry(i) = (rysc(i) +drysc(i))/sc(i)/rr(i)
         by(i) = (bysc(i) +dbysc(i))/sc(i)*rr(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
         vy(i) = ry(i)/ro(i)
         vv=vx(i)**2+vy(i)**2
c         bb=by(i)**2
         bb=bx(i)**2+by(i)**2
         pr(i) = (gm-1)*(ee(i) - 0.5*ro(i)*vv - pi8i*bb)
      enddo

      return
      end
