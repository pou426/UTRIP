c======================================================================|
      subroutine mlw_sm(ro,pr,vx,vy,by,bx,bxm,dt,qav,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_sm
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * special relativistic 
c
c INPUTS & OUTPUTS
c
c OUTPUTS
c    None
c
c INPUTS
c
c HISTORY
c    written 2004-2-18 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)

      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
      dimension gl(ix),de(ix),ee(ix),rx(ix),ry(ix)
      dimension ez(ix)
      dimension pm(ix),el(ix)
      dimension pmh(ix),elh(ix)

      dimension bx(ix),bxm(ix)

      dimension roh(ix),eeh(ix)
      dimension prh(ix),vxh(ix),vyh(ix),glh(ix)
      dimension deh(ix),rxh(ix),ryh(ix),byh(ix)
      dimension ezh(ix)

      dimension dde(ix),dee(ix),drx(ix),dry(ix),dby(ix)

      dimension fx(ix),qx(ix)

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
c----------------------------------------------------------------------|
c     initialize dde etc.                                   
c----------------------------------------------------------------------|
      
      do i=1,ix
         dde(i) = 0.0
         dee(i) = 0.0
         drx(i) = 0.0
         dry(i) = 0.0
         dby(i) = 0.0
         glh(i) = 0.d0
         roh(i) = 0.d0
         prh(i) = 0.d0
         vxh(i) = 0.d0
         vyh(i) = 0.d0
         byh(i) = 0.d0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix       
         gl(i) = 1./sqrt(1-vx(i)**2-vy(i)**2)
         ez(i) = -vx(i)*by(i)+vy(i)*bx(i)
         de(i) = gl(i)*ro(i)
         ef2 = ez(i)**2
         bf2 = bx(i)**2+by(i)**2
         pm(i) = pi8i*(bf2+ef2)
         en = ro(i)+pr(i)/(gm-1)
         el(i) = (en+pr(i))*gl(i)**2
         ee(i) = el(i)-pr(i)-de(i)+pm(i)
         rx(i) = el(i)*vx(i)+pi4i*(-ez(i)*by(i))
         ry(i) = el(i)*vy(i)+pi4i*( ez(i)*bx(i))
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix       
         fx(i)= de(i)*vx(i)
      enddo
      call mlwhalf(de,deh,dde,dt,fx,dxi,dxim,ix)

c---  energy ---
      do i=1,ix
         fx(i)= rx(i)-de(i)*vx(i)
      enddo
      call mlwhalf(ee,eeh,dee,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix
         fx(i)= pr(i)+pm(i)+el(i)*vx(i)**2-pi4i*(bx(i)**2)
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix)

c---  y-momentum ---
      do i=1,ix
         fx(i)= el(i)*vx(i)*vy(i)-pi4i*(bx(i)*by(i))
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix)

c---  y-magnetic ---
      do i=1,ix
         fx(i)= -ez(i)
      enddo
      call mlwhalf(by,byh,dby,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      mix=100
      tolf=1.d-20
      tolx=1.d-20
c---- initial guess
      do i=1,ix
        roh(i)=ro(i)
        prh(i)=pr(i)
        vxh(i)=vx(i)
        vyh(i)=vy(i)
      enddo
c---- solve nonlinear equations
      call rtnewt_sm(roh,prh,vxh,vyh,glh
     &     ,deh,eeh,rxh,ryh,byh,bxm,gm,mix,tolf,tolx,ix,1,ix-1)

      do i=1,ix-1
        ezh(i)=-vxh(i)*byh(i)+vyh(i)*bxm(i)
        ef2 = ezh(i)**2
        bf2 = bxm(i)**2+byh(i)**2
        pmh(i)= pi8i*(bf2+ef2) 
        enh = roh(i)+prh(i)/(gm-1)
        elh(i) = (enh+prh(i))*glh(i)**2
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= deh(i)*vxh(i)
      enddo
      call mlwfull(dde ,dt,fx,dxi,ix)

c---  energy     ---
      do i=1,ix-1
         fx(i)= rxh(i)-deh(i)*vxh(i)
      enddo
      call mlwfull(dee ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         fx(i)= prh(i)+pmh(i)+elh(i)*vxh(i)**2-pi4i*(bxm(i)**2)
      enddo
      call mlwfull(drx,dt,fx,dxi,ix)

c---  y-momentum ---
      do i=1,ix-1
         fx(i)= elh(i)*vxh(i)*vyh(i)-pi4i*(bxm(i)*byh(i))
      enddo
      call mlwfull(dry,dt,fx,dxi,ix)

c---  y-magnetic ---
      do i=1,ix-1
         fx(i)= -ezh(i)
      enddo
      call mlwfull(dby,dt,fx,dxi,ix)
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
c     qav=3.0d0
      zero=0.0d0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
      call mlwartv(de,dde,dt,qx,dxi,dxim,ix)
      call mlwartv(ee,dee,dt,qx,dxi,dxim,ix)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix)
      call mlwartv(by,dby,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         de(i) = de(i) +dde(i)
         ee(i) = ee(i) +dee(i)
         rx(i) = rx(i) +drx(i)
         ry(i) = ry(i) +dry(i)
         by(i) = by(i) +dby(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|

c---- solve nonlinear equations
      call rtnewt_sm(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,by,bx,gm,mix,tolf,tolx,ix,2,ix-1)


      return
      end
