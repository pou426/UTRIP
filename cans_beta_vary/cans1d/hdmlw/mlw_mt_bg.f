c======================================================================|
      subroutine mlw_mt_bg(ro,vx,vy,by,bx,bxm,dt,qav,cs2
     &             ,gx,gxm,sc,dsc,scm,dscm,rr,rrm,drr,drrm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_mt_bg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal MHD
c        * axial symmetry
c        * non-uniform poloidal magnetic field
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
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

      dimension ro(ix),vx(ix),vy(ix),by(ix)
      dimension rx(ix),ry(ix),ez(ix)

      dimension roh(ix),vxh(ix),vyh(ix),byh(ix)
      dimension ezh(ix)

      dimension bx(ix),bxm(ix)

      dimension sc(ix),scm(ix),dsc(ix),dscm(ix)
      dimension rosc(ix),rxsc(ix),rysc(ix),bysc(ix)
      dimension rosch(ix),rxsch(ix),rysch(ix),bysch(ix)
      dimension drosc(ix),drxsc(ix),drysc(ix),dbysc(ix)

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
         drxsc(i) = 0.0
         drysc(i) = 0.0
         dbysc(i) = 0.0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         rx(i) = ro(i)*vx(i)
         ry(i) = ro(i)*vy(i)
      enddo
      do i=1,ix
         ez(i) = -vx(i)*by(i)+vy(i)*bx(i)
      enddo
      do i=1,ix       
         rosc(i) = ro(i)*sc(i)
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

c---  x-momentum ---
      do i=1,ix
         bbm=by(i)**2
c        bbm=by(i)**2-bx(i)**2
         bb =by(i)**2
c        bb =by(i)**2+bx(i)**2
         fx(i)= ro(i)*vx(i)**2+cs2*ro(i)+bbm*pi8i
         ss(i)= ro(i)*gx(i)
     &         +(ro(i)*vy(i)**2-by(i)**2*pi4i)/rr(i)*drr(i)
         fx(i)= fx(i)*sc(i)
         ss(i)= ss(i)*sc(i)+(cs2*ro(i)+bb*pi8i)*dsc(i)
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
         vxh(i)   = rxsch(i)/rosch(i)
         vyh(i)   = rysch(i)/rosch(i)/rrm(i)
         byh(i)   = bysch(i)/scm(i)*rrm(i)
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

c---  x-momentum ---
      do i=1,ix-1
         bbm=byh(i)**2
         bb =byh(i)**2
c        bbm=byh(i)**2-bxm(i)**2
c        bb =byh(i)**2+bxm(i)**2
         fx(i)= roh(i)*vxh(i)**2+cs2*roh(i)+bbm*pi8i
         ss(i)= roh(i)*gxm(i)
     &        +(roh(i)*vyh(i)**2-byh(i)**2*pi4i)/rrm(i)*drrm(i)
         fx(i)= fx(i)*scm(i)
         ss(i)= ss(i)*scm(i)+(cs2*roh(i)+bb*pi8i)*dscm(i)
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
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix)
      call mlwartv(rysc,drysc,dt,qx,dxi,dxim,ix)
      call mlwartv(bysc,dbysc,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = (rosc(i) +drosc(i))/sc(i)
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
      enddo

      return
      end
