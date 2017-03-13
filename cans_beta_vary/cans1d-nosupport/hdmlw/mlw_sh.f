c======================================================================|
      subroutine mlw_sh(ro,pr,vx,dt,qav,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_r
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * general relativistic hydrodynamics
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity along the x-cordinate
c    gl(ix), glh(ix): [double] Lorentz factor
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ix: [integer] dimension size
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    hh(ix,4), hhm(ix,4) : [double] metric, h0,h1,h2,h3
c    gg(ix,3), ggm(ix,3) : [double] curvature, G01, G21, G31
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)

      dimension ro(ix),pr(ix),vx(ix)
      dimension gl(ix),glh(ix),de(ix),ee(ix),rx(ix)

      dimension roh(ix),eeh(ix)
      dimension prh(ix),vxh(ix)
      dimension deh(ix),rxh(ix)

      dimension dde(ix),dee(ix),drx(ix)

      dimension fx(ix),qx(ix)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.d0/dxm(i)
         dxi(i) = 1.d0/dx(i)
      enddo
c----------------------------------------------------------------------|
c     initialize dde etc.                                   
c----------------------------------------------------------------------|
      do i=1,ix
         dde(i) = 0.d0
         drx(i) = 0.d0
         dee(i) = 0.d0
         glh(i) = 0.d0
         roh(i) = 0.d0
         prh(i) = 0.d0
         vxh(i) = 0.d0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         gl(i) = 1.d0/sqrt(1.d0-vx(i)**2)
         de(i) = gl(i)*ro(i)
         en = ro(i)+pr(i)/(gm-1.d0)
         ee(i) = (en+pr(i))*gl(i)**2-pr(i)-de(i)
         rx(i) = (en+pr(i))*gl(i)**2*vx(i)
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
         en = ro(i)+pr(i)/(gm-1)
         fx(i)= pr(i)+(en+pr(i))*gl(i)**2*vx(i)**2
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      tolx=1.d-10
      tolf=1.d-10
      mix=1000
      glmax=1.d1

c     call rtbis_sh(roh,prh,vxh,glh,deh,eeh,rxh,gm,mix,tolf,tolx,glmax
c    &       ,ix,1,ix-1)

      call rtnewt_sh(roh,prh,vxh,glh,deh,eeh,rxh,gm,mix,tolf,tolx
     &       ,ix,1,ix-1)

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
         enh = roh(i)+prh(i)/(gm-1)
         fx(i)= prh(i)+(enh+prh(i))*glh(i)**2*vxh(i)**2
      enddo
      call mlwfull(drx,dt,fx,dxi,ix)
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
c     qav=4.0d0
      zero=0.0d0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
      call mlwartv(de,dde,dt,qx,dxi,dxim,ix)
      call mlwartv(ee,dee,dt,qx,dxi,dxim,ix)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         de(i) = de(i) +dde(i)
         rx(i) = rx(i) +drx(i)
         ee(i) = ee(i) +dee(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|

c     call rtbis_sh(ro,pr,vx,gl,de,ee,rx,gm,mix,tolf,tolx,glmax
c    &       ,ix,2,ix-1)

      call rtnewt_sh(ro,pr,vx,gl,de,ee,rx,gm,mix,tolf,tolx
     &       ,ix,2,ix-1)

      return
      end
