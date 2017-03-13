c======================================================================|
      subroutine mlw_sh(ro,pr,vx,vy,dt,qav,gm,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  mlw_sh
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * special relativity
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    az(ix,jx): [double] magnetic vector potential
c 
c OUTPUTS
c    None
c 
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c    
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dt: [double] delta time 
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx)
      dimension dde(ix,jx),dee(ix,jx),drx(ix,jx),dry(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)

      dimension gl(ix,jx),de(ix,jx),el(ix,jx)
      dimension glm(ix,jx),deh(ix,jx),elh(ix,jx)
c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxi(i) = 1.0/dx(i)
         dxim(i) = 1.0/dxm(i)
      enddo
      do i=2,ix-1
         ux1(i)  = 0.5*dxm(i-1)/dx(i)
         ux0(i)  = 0.5*dxm(i)/dx(i)
      enddo

      do j=1,jx
         dyi(j)  = 1.0/dy(j)
         dyim(j) = 1.0/dym(j)
      enddo
      do j=2,jx-1
         uy1(j)  = 0.5*dym(j-1)/dy(j)
         uy0(j)  = 0.5*dym(j)/dy(j)
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         dde(i,j)= 0.d0
         dee(i,j)= 0.d0
         drx(i,j)= 0.d0
         dry(i,j)= 0.d0
         glm(i,j)= 0.d0
         roh(i,j)= 0.d0
         prh(i,j)= 0.d0
         vxh(i,j)= 0.d0
         vyh(i,j)= 0.d0
      enddo
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         gl(i,j) = 1./sqrt(1-vx(i,j)**2-vy(i,j)**2)
         de(i,j) = gl(i,j)*ro(i,j)
         en = ro(i,j)+pr(i,j)/(gm-1)
         el(i,j) = (en+pr(i,j))*gl(i,j)**2
         ee(i,j) = el(i,j)-pr(i,j)-de(i,j)
         rx(i,j) = el(i,j)*vx(i,j)
         ry(i,j) = el(i,j)*vy(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= de(i,j)*vx(i,j)
         fy(i,j)= de(i,j)*vy(i,j)
      enddo
      enddo
      call mlwhalf(de ,deh ,dde,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  energy ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= rx(i,j)-de(i,j)*vx(i,j)
         fy(i,j)= ry(i,j)-de(i,j)*vy(i,j)
      enddo
      enddo
      call mlwhalf(ee ,eeh ,dee ,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  x-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= el(i,j)*vx(i,j)**2+pr(i,j)
         fy(i,j)= el(i,j)*vx(i,j)*vy(i,j)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= el(i,j)*vx(i,j)*vy(i,j)
         fy(i,j)= el(i,j)*vy(i,j)**2+pr(i,j)
      enddo
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      mix=100
      tolf=1.d-20
      tolx=1.d-20
c---- initial guess
      do j=1,jx
      do i=1,ix
        roh(i,j)=ro(i,j)
        prh(i,j)=pr(i,j)
        vxh(i,j)=vx(i,j)
        vyh(i,j)=vy(i,j)
      enddo
      enddo
c---- solve nonlinear equations
      call rtnewt_sh(roh,prh,vxh,vyh,glm
     &     ,deh,eeh,rxh,ryh,gm,mix,tolf,tolx
     &     ,ix,1,ix-1,jx,1,jx-1)

      do j=1,jx-1
      do i=1,ix-1
        enh = roh(i,j)+prh(i,j)/(gm-1)
        elh(i,j) = (enh+prh(i,j))*glm(i,j)**2
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= deh(i,j)*vxh(i,j)
         fy(i,j)= deh(i,j)*vyh(i,j)
      enddo
      enddo
      call mlwfull(dde ,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  energy     ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= rxh(i,j)-deh(i,j)*vxh(i,j)
         fy(i,j)= ryh(i,j)-deh(i,j)*vyh(i,j)
      enddo
      enddo
      call mlwfull(dee ,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= elh(i,j)*vxh(i,j)**2+prh(i,j)
         fy(i,j)= elh(i,j)*vxh(i,j)*vyh(i,j)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= elh(i,j)*vxh(i,j)*vyh(i,j)
         fy(i,j)= elh(i,j)*vyh(i,j)**2+prh(i,j)
      enddo
      enddo
      call mlwfull(dry,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c-------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0
      do j=1,jx-1
      do i=1,ix-1
         qx(i,j)=qav*dxm(i)*max(zero,abs(vx(i+1,j)-vx(i,j))-1.0e-4)
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
         qy(i,j)=qav*dym(j)*max(zero,abs(vy(i,j+1)-vy(i,j))-1.0e-4)
      enddo
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(de,dde,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(ee,dee,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         de(i,j) = de(i,j) +dde(i,j)
         ee(i,j) = ee(i,j) +dee(i,j)
         rx(i,j) = rx(i,j) +drx(i,j)
         ry(i,j) = ry(i,j) +dry(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
c---- solve nonlinear equations
      call rtnewt_sh(ro,pr,vx,vy,gl
     &     ,de,ee,rx,ry,gm,mix,tolf,tolx
     &     ,ix,2,ix-1,jx,2,jx-1)


      return
      end
