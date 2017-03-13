c======================================================================|
      subroutine mlw_sh_o(ro,pr,vx,dt,qav,gm,hh,hhm,gg,ggm,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_sh_o
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * special relativistic hydrodynamics
c        * orthogonal general coordinate
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity along the x-cordinate
c    gl(ix), glm(ix): [double] Lorentz factor
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
c    hh(ix,3), hhm(ix,3) : [double] metric, h0,h1,h2,h3
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)

      dimension ro(ix),pr(ix),vx(ix)
      dimension gl(ix),glm(ix),de(ix),ee(ix),rx(ix),el(ix)

      dimension roh(ix),eeh(ix),prh(ix),vxh(ix),deh(ix),rxh(ix),elh(ix)

      dimension hh(ix,3),hhm(ix,3)
      dimension gg(ix,3),ggm(ix,3)
      dimension desc(ix),eesc(ix),rxsc(ix)
      dimension desch(ix),eesch(ix),rxsch(ix)
      dimension ddesc(ix),deesc(ix),drxsc(ix)

      dimension fx(ix),qx(ix),ss(ix)
      dimension sc(ix)

c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.d0/dxm(i)
         dxi(i) = 1.d0/dx(i)
      enddo
      do i=2,ix-1
         ux1(i)  = 0.5*dxm(i-1)/dx(i)
         ux0(i)  = 0.5*dxm(i)/dx(i)
      enddo
c----------------------------------------------------------------------|
c     initialize dde etc.                                   
c----------------------------------------------------------------------|
      do i=1,ix
         ddesc(i) = 0.d0
         drxsc(i) = 0.d0
         deesc(i) = 0.d0
         glm(i) = 0.d0
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
         el(i) = (en+pr(i))*gl(i)**2
         ee(i) = el(i)-pr(i)-de(i)
         rx(i) = el(i)*vx(i)
      enddo
      do i=1,ix       
         sc(i)=hh(i,1)*hh(i,2)*hh(i,3)
         desc(i) = de(i)*sc(i)
         rxsc(i) = rx(i)*sc(i)
         eesc(i) = ee(i)*sc(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix       
         fx(i)= de(i)*vx(i)
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
      enddo
      call mlwhalf(desc,desch,ddesc,dt,fx,dxi,dxim,ix)

c---  energy ---
      do i=1,ix
         fx(i)= rx(i)-de(i)*vx(i)
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
      enddo
      call mlwhalf(eesc,eesch,deesc,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix
         t11= pr(i)+el(i)*vx(i)**2
         t22= pr(i)
         t33= pr(i)
         t12= el(i)*vx(i)*vy(i)
         t23= el(i)*vy(i)*vz(i)
         t31= el(i)*vz(i)*vx(i)
         fx(i)= t11
         fx(i)= fx(i)*hh(i,2)*hh(i,3)
         ss(i)= gg(i,3)*t22+gg(i,2)*t33-gg(i,3)*t12-gg(i,2)*t31
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix)
      call mlwsrch(rxsch,drxsc,dt,ss,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         scm=hhm(i,1)*hhm(i,2)*hhm(i,3)
         deh(i) = desch(i)/scm
         eeh(i) = eesch(i)/scm
         rxh(i) = rxsch(i)/scm
      enddo

      tolx=1.d-10
      tolf=1.d-10
      mix=1000
      glmax=1.d1

c     call rtbis_sh(roh,prh,vxh,glm,deh,eeh,rxh,gm,mix,tolf,tolx,glmax
c    &       ,ix,1,ix-1)

      call rtnewt_sh(roh,prh,vxh,glm,deh,eeh,rxh,gm,mix,tolf,tolx
     &       ,ix,1,ix-1)

      do i=1,ix-1
        enh = roh(i)+prh(i)/(gm-1)
        elh(i) = (enh+prh(i))*glm(i)**2
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= deh(i)*vxh(i)
         fx(i)= fx(i)*hhm(i,2)*hhm(i,3)
      enddo
      call mlwfull(ddesc ,dt,fx,dxi,ix)

c---  energy     ---
      do i=1,ix-1
         fx(i)= rxh(i)-deh(i)*vxh(i)
      enddo
      call mlwfull(deesc ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         t11= prh(i)+elh(i)*vxh(i)**2
         t22= prh(i)
         t33= prh(i)
         t12= elh(i)*vxh(i)*vyh(i)
         t23= elh(i)*vyh(i)*vzh(i)
         t31= elh(i)*vzh(i)*vxh(i)
         fx(i)= t11
         fx(i)= fx(i)*hhm(i,2)*hhm(i,3)
         ss(i)= ggm(i,3)*t22+ggm(i,2)*t33-ggm(i,3)*t12-ggm(i,2)*t31
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ix)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix)
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
c     qav=3.0d0
      zero=0.0d0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
      call mlwartv(desc,ddesc,dt,qx,dxi,dxim,ix)
      call mlwartv(eesc,deesc,dt,qx,dxi,dxim,ix)
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         de(i) = (desc(i) +ddesc(i))/sc(i)
         rx(i) = (rxsc(i) +drxsc(i))/sc(i)
         ee(i) = (eesc(i) +deesc(i))/sc(i)
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
