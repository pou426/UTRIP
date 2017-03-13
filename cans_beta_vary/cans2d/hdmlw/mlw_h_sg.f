c======================================================================|
      subroutine mlw_h_sg(ro,pr,vx,vy,dt,qav,gm
     &                  ,gx,gxm,gy,gym,x,xm,dx,dxm,ix,y,ym,dy,dym,jx)
c======================================================================|
c
c NAME  mlw_h_sg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * hydrodynamics
c        * Spherical coordinate, axis-symmetry
c        * gravity
c 
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c 
c    gx(ix,jx), gxm(ix,jx) : [double] gravity
c    gy(ix,jx), gym(ix,jx) : [double] gravity
c    x(ix), xm(ix): [double] coordinate
c    y(jx), ym(jx): [double] coordinate
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
      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx)
      dimension dyi(jx),dyim(jx)
      dimension uy0(jx),uy1(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension prh(ix,jx),roh(ix,jx),eeh(ix,jx),vxh(ix,jx),vyh(ix,jx)
      dimension rosc(ix,jx),eesc(ix,jx),rxsc(ix,jx),rysc(ix,jx)
      dimension rosch(ix,jx),eesch(ix,jx),rxsch(ix,jx),rysch(ix,jx)
      dimension drosc(ix,jx),deesc(ix,jx),drxsc(ix,jx),drysc(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)
      dimension ss(ix,jx)
      dimension x(ix),xm(ix)
      dimension y(jx),ym(jx)
      dimension siny(jx),sinym(jx)
      dimension gx(ix,jx), gxm(ix,jx)
      dimension gy(ix,jx), gym(ix,jx)
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

      do j=1,jx
         siny(j) = sin(y(j))
         sinym(j) = sin(ym(j))
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         drosc(i,j) = 0.0
         deesc(i,j) = 0.0
         drxsc(i,j) = 0.0
         drysc(i,j) = 0.0
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vy(i,j)**2
         ee(i,j) = pr(i,j)/(gm-1)+0.5*ro(i,j)*vv
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
      enddo
      enddo
      do j=1,jx
      do i=1,ix
         rosc(i,j) = ro(i,j)*x(i)**2*siny(j)
         eesc(i,j) = ee(i,j)*x(i)**2*siny(j)
         rxsc(i,j) = rx(i,j)*x(i)**2*siny(j)
         rysc(i,j) = ry(i,j)*x(i)*siny(j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)
         fx(i,j)= fx(i,j)*x(i)**2*siny(j)
         fy(i,j)= fy(i,j)*x(i)*siny(j)
      enddo
      enddo
      call mlwhalf(rosc,rosch,drosc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  energy ---
      do j=1,jx
      do i=1,ix
c        vv=vx(i,j)**2+vy(i,j)**2
c        ep=pr(i,j)*gm/(gm-1.)+0.5*ro(i,j)*vv
         ep=pr(i,j)+ee(i,j)
         fx(i,j)= ep*vx(i,j)
         fy(i,j)= ep*vy(i,j)
         ss(i,j)= ro(i,j)*(vx(i,j)*gx(i,j)+vy(i,j)*gy(i,j))
         fx(i,j)= fx(i,j)*x(i)**2*siny(j)
         fy(i,j)= fy(i,j)*x(i)*siny(j)
         ss(i,j)= ss(i,j)*x(i)**2*siny(j)
      enddo
      enddo
      call mlwhalf(eesc,eesch,deesc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(eesch,deesc,dt,ss,ix,jx)

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)
         fy(i,j)= ro(i,j)*vx(i,j)*vy(i,j)
         ss(i,j)= ro(i,j)*gx(i,j)
         fx(i,j)= fx(i,j)*x(i)**2*siny(j)
         fy(i,j)= fy(i,j)*x(i)*siny(j)
         ss(i,j)= ss(i,j)*x(i)**2*siny(j)
     &            +(ro(i,j)*vy(i,j)**2+2*pr(i,j))*x(i)*siny(j)
      enddo
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rxsch,drxsc,dt,ss,ix,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)**2+pr(i,j)
         ss(i,j)= ro(i,j)*gy(i,j)
         fx(i,j)= fx(i,j)*x(i)*siny(j)
         fy(i,j)= fy(i,j)*siny(j)
         ss(i,j)= ss(i,j)*x(i)*siny(j)
     &           +(ro(i,j)*vx(i,j)*vy(i,j))*(-2)*siny(j)
     &          + pr(i,j)*cos(y(j))
      enddo
      enddo
      call mlwhalf(rysc,rysch,drysc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rysch,drysc,dt,ss,ix,jx)
c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         roh(i,j)   = rosch(i,j)/xm(i)**2/sinym(j)
         eeh(i,j)   = eesch(i,j)/xm(i)**2/sinym(j)
         vxh(i,j)   = rxsch(i,j)/roh(i,j)/xm(i)**2/sinym(j)
         vyh(i,j)   = rysch(i,j)/roh(i,j)/xm(i)/sinym(j)
         vv=vxh(i,j)**2+vyh(i,j)**2
         prh(i,j)   = (gm-1)*(eeh(i,j)-0.5*roh(i,j)*vv)
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)
         fy(i,j)= roh(i,j)*vyh(i,j)
         fx(i,j)= fx(i,j)*xm(i)**2*sinym(j)
         fy(i,j)= fy(i,j)*xm(i)*sinym(j)
      enddo
      enddo
      call mlwfull(drosc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  energy     ---
      do j=1,jx-1
      do i=1,ix-1
c        vv=vxh(i,j)**2+vyh(i,j)**2
c        eph=prsch(i,j)*gm/(gm-1.)+0.5*rosch(i,j)*vv
         eph=prh(i,j)+eeh(i,j)
         fx(i,j)= eph*vxh(i,j)
         fy(i,j)= eph*vyh(i,j)
         ss(i,j)= roh(i,j)*(vxh(i,j)*gxm(i,j)+vyh(i,j)*gym(i,j))
         fx(i,j)= fx(i,j)*xm(i)**2*sinym(j)
         fy(i,j)= fy(i,j)*xm(i)*sinym(j)
         ss(i,j)= ss(i,j)*xm(i)**2*sinym(j)
      enddo
      enddo
      call mlwfull(deesc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(deesc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+prh(i,j)
         fy(i,j)= roh(i,j)*vxh(i,j)*vyh(i,j)
         ss(i,j)= roh(i,j)*gxm(i,j)
         fx(i,j)= fx(i,j)*xm(i)**2*sinym(j)
         fy(i,j)= fy(i,j)*xm(i)*sinym(j)
         ss(i,j)= ss(i,j)*xm(i)**2*sinym(j)
     &            +(roh(i,j)*vyh(i,j)**2+2*prh(i,j))*xm(i)*sinym(j)
      enddo
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vyh(i,j)*vxh(i,j)
         fy(i,j)= roh(i,j)*vyh(i,j)**2+prh(i,j)
         ss(i,j)= roh(i,j)*gym(i,j)
         fx(i,j)= fx(i,j)*xm(i)*sinym(j)
         fy(i,j)= fy(i,j)*sinym(j)
         ss(i,j)= ss(i,j)*xm(i)*sinym(j)
     &           +(roh(i,j)*vxh(i,j)*vyh(i,j))*(-2)*sinym(j)
     &          + prh(i,j)*cos(ym(j))
      enddo
      enddo
      call mlwfull(drysc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drysc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)
c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=1.0
      zero=0.0e0
      do j=1,jx
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
      call mlwartv(rosc,drosc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(eesc,deesc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rysc,drysc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         rosc(i,j) = rosc(i,j) +drosc(i,j)
         eesc(i,j) = eesc(i,j) +deesc(i,j)
         rxsc(i,j) = rxsc(i,j) +drxsc(i,j)
         rysc(i,j) = rysc(i,j) +drysc(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=2,jx-1  
      do i=2,ix-1  
         ro(i,j) = rosc(i,j)/x(i)**2/siny(j)
         ee(i,j) = eesc(i,j)/x(i)**2/siny(j)
         vx(i,j) = rxsc(i,j)/ro(i,j)/x(i)**2/siny(j)
         vy(i,j) = rysc(i,j)/ro(i,j)/x(i)/siny(j)
         vv=vx(i,j)**2+vy(i,j)**2
         pr(i,j) = (gm-1)*(ee(i,j) - 0.5*ro(i,j)*vv)
      enddo
      enddo

      return
      end
