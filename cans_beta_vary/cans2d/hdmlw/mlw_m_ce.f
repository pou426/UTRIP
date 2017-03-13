c======================================================================|
      subroutine mlw_m_ce(ro,pr,vx,vz,bx,bz,ay
     &                 ,dt,qav,gm,et,etm,x,xm,dx,dxm,ix,dz,dzm,jx)
c======================================================================|
c
c NAME  mlw_m_ce
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * Cylindrical coordinate, axis-symmetry
c        * resistivity
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
c    et(ix), etm(ix): [double] resistivity
c    dx(ix), dxm(ix): [double] grid spacing
c    dz(jx), dzm(jx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c 
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ux0(ix),ux1(ix)
      dimension x(ix),xm(ix)
      dimension dz(jx),dzm(jx)
      dimension dzi(jx),dzim(jx)
      dimension uz0(jx),uz1(jx)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vz(ix,jx)
     &         ,bx(ix,jx),bz(ix,jx)
      dimension ey(ix,jx)
      dimension ay(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),rz(ix,jx)

      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),rzh(ix,jx)
     &         ,bxh(ix,jx),bzh(ix,jx)
      dimension eyh(ix,jx)
      dimension ayh(ix,jx)

      dimension prh(ix,jx),vxh(ix,jx),vzh(ix,jx)
      dimension dro(ix,jx),dee(ix,jx),drx(ix,jx),drz(ix,jx)
     &         ,dbx(ix,jx),dbz(ix,jx)
     &         ,day(ix,jx)

      dimension fx(ix,jx),qx(ix,jx)
      dimension fz(ix,jx),qz(ix,jx)
      dimension ss(ix,jx)

      dimension cy(ix,jx), cyh(ix,jx)
      dimension et(ix,jx), etm(ix,jx)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
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
         dzi(j)  = 1.0/dz(j)
         dzim(j) = 1.0/dzm(j)
      enddo
      do j=2,jx-1
         uz1(j)  = 0.5*dzm(j-1)/dz(j)
         uz0(j)  = 0.5*dzm(j)/dz(j)
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         dro(i,j)= 0.0
         dee(i,j)= 0.0
         drx(i,j)= 0.0
         drz(i,j)= 0.0
         dbx(i,j)= 0.0
         dbz(i,j)= 0.0
         day(i,j)= 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vz(i,j)**2
         bb=bx(i,j)**2+bz(i,j)**2
         ee(i,j) = pr(i,j)/(gm-1)+0.5*ro(i,j)*vv+pi8i*bb
         rx(i,j) = ro(i,j)*vx(i,j)
         rz(i,j) = ro(i,j)*vz(i,j)
      enddo
      enddo

      call bbtocy_c(cy,bz,bx,dz,dx,ix,jx)

      do j=1,jx
      do i=1,ix
         ey(i,j) = -vz(i,j)*bx(i,j)+vx(i,j)*bz(i,j)+et(i,j)*cy(i,j)
      enddo
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)
         fz(i,j)= ro(i,j)*vz(i,j)
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwhalf(ro ,roh ,dro,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(roh ,dro ,dt,ss,ix,jx)

c---  energy ---
      do j=1,jx
      do i=1,ix       
         vv=vx(i,j)**2+vz(i,j)**2
         ep    = pr(i,j)*gm/(gm-1.)+0.5*ro(i,j)*vv
         fx(i,j)= ep*vx(i,j) +(bz(i,j)*ey(i,j))*pi4i
         fz(i,j)= ep*vz(i,j) +(-bx(i,j)*ey(i,j))*pi4i
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwhalf(ee ,eeh ,dee ,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(eeh ,dee ,dt,ss,ix,jx)

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)
     &          +pi8i*(bz(i,j)**2-bx(i,j)**2)
         fz(i,j)= ro(i,j)*vx(i,j)*vz(i,j)-pi4i*bx(i,j)*bz(i,j)
         ss(i,j)= -(ro(i,j)*(vx(i,j)**2)
     &              +pi4i*(-bx(i,j)**2))/x(i)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rxh ,drx ,dt,ss,ix,jx)

c---  z-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= ro(i,j)*vz(i,j)*vx(i,j)-pi4i*bz(i,j)*bx(i,j)
         fz(i,j)= ro(i,j)*vz(i,j)**2+pr(i,j)
     &          +pi8i*(bx(i,j)**2-bz(i,j)**2)
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwhalf(rz,rzh,drz,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rzh ,drz ,dt,ss,ix,jx)

c---  x-magnetic ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= 0.
         fz(i,j)= -ey(i,j)
      enddo
      enddo
      call mlwhalf(bx,bxh,dbx,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)

c---  z-magnetic ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= ey(i,j)
         fz(i,j)= 0.
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwhalf(bz,bzh,dbz,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(bzh ,dbz ,dt,ss,ix,jx)

c---  y-magnetic potential ---
      do j=1,jx
      do i=1,ix
         ss(i,j)= -ey(i,j)
      enddo
      enddo
      call mlwsrch(ayh ,day ,dt,ss,ix,jx)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         vxh(i,j)   = rxh(i,j)/roh(i,j)
         vzh(i,j)   = rzh(i,j)/roh(i,j)
         vv=vxh(i,j)**2+vzh(i,j)**2
         bb=bxh(i,j)**2+bzh(i,j)**2
         prh(i,j)   = (gm-1)*(eeh(i,j)-0.5*roh(i,j)*vv-pi8i*bb)
      enddo
      enddo

      call bbtocy_c(cyh,bzh,bxh,dzm,dxm,ix,jx)

      do j=1,jx
      do i=1,ix
        eyh(i,j)=-vzh(i,j)*bxh(i,j)+vxh(i,j)*bzh(i,j)+etm(i,j)*cyh(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)
         fz(i,j)= roh(i,j)*vzh(i,j)
         ss(i,j)= -fx(i,j)/xm(i)
      enddo
      enddo
      call mlwfull(dro ,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(dro,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  energy     ---
      do j=1,jx-1
      do i=1,ix-1
         vv=vxh(i,j)**2+vzh(i,j)**2
         ep    = prh(i,j)*gm/(gm-1.)+0.5*roh(i,j)*vv
         fx(i,j)= ep*vxh(i,j)
     &            +( bzh(i,j)*eyh(i,j))*pi4i
         fz(i,j)= ep*vzh(i,j)
     &            +(-bxh(i,j)*eyh(i,j))*pi4i
         ss(i,j)= -fx(i,j)/xm(i)
      enddo
      enddo
      call mlwfull(dee ,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(dee,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+prh(i,j)
     &          +pi8i*(bzh(i,j)**2-bxh(i,j)**2)
         fz(i,j)= roh(i,j)*vxh(i,j)*vzh(i,j)-pi4i*bxh(i,j)*bzh(i,j)
         ss(i,j)= -(roh(i,j)*(vxh(i,j)**2)
     &              +pi4i*(-bxh(i,j)**2))/xm(i)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drx,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  z-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vzh(i,j)*vxh(i,j)-pi4i*bzh(i,j)*bxh(i,j)
         fz(i,j)= roh(i,j)*vzh(i,j)**2+prh(i,j)
     &          +pi8i*(bxh(i,j)**2-bzh(i,j)**2)
         ss(i,j)= -fx(i,j)/xm(i)
      enddo
      enddo
      call mlwfull(drz,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drz,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  x-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= 0.
         fz(i,j)= -eyh(i,j)
      enddo
      enddo
      call mlwfull(dbx,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)

c---  z-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= eyh(i,j)
         fz(i,j)= 0.
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwfull(dbz,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(dbz,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  y-magnetic potential ---
      do j=1,jx
      do i=1,ix
         ss(i,j)= -eyh(i,j)
      enddo
      enddo
      call mlwsrcf(day,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

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
         qz(i,j)=qav*dzm(j)*max(zero,abs(vz(i,j+1)-vz(i,j))-1.0e-4)
      enddo
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(ro,dro,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(ee,dee,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(rz,drz,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(bx,dbx,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(bz,dbz,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         ro(i,j) = ro(i,j) +dro(i,j)
         ee(i,j) = ee(i,j) +dee(i,j)
         rx(i,j) = rx(i,j) +drx(i,j)
         rz(i,j) = rz(i,j) +drz(i,j)
         bx(i,j) = bx(i,j) +dbx(i,j)
         bz(i,j) = bz(i,j) +dbz(i,j)
         ay(i,j) = ay(i,j) +day(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1  
         vx(i,j) = rx(i,j)/ro(i,j)
         vz(i,j) = rz(i,j)/ro(i,j)
         vv=vx(i,j)**2+vz(i,j)**2
         bb=bx(i,j)**2+bz(i,j)**2
         pr(i,j) = (gm-1)*(ee(i,j) - 0.5*ro(i,j)*vv - pi8i*bb)
      enddo
      enddo

      return
      end
