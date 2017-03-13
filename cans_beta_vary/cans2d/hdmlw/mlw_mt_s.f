c======================================================================|
      subroutine mlw_mt_s(ro,vx,vy,bx,by,az,dt,qav,cs2
     &                  ,x,xm,dx,dxm,ix,y,ym,dy,dym,jx)
c======================================================================|
c
c NAME  mlw_mt_s
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal MHD
c        * Spherical coordinate, axis-symmetry
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
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
c    x(ix), xm(ix): [double] coordinate
c    y(jx), ym(jx): [double] coordinate
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dt: [double] delta time
c    cs2: [double] square of sound speed
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
     &         ,bx(ix,jx),by(ix,jx)
      dimension rx(ix,jx),ry(ix,jx),ez(ix,jx)
      dimension az(ix,jx)
      dimension prh(ix,jx),roh(ix,jx),vxh(ix,jx),vyh(ix,jx)
     &         ,bxh(ix,jx),byh(ix,jx),ezh(ix,jx)
     &         ,azh(ix,jx)

      dimension rosc(ix,jx),rxsc(ix,jx),rysc(ix,jx)
     &         ,bxsc(ix,jx),bysc(ix,jx)
      dimension rosch(ix,jx),rxsch(ix,jx),rysch(ix,jx)
     &         ,bxsch(ix,jx),bysch(ix,jx)
      dimension drosc(ix,jx),drxsc(ix,jx),drysc(ix,jx)
     &         ,dbxsc(ix,jx),dbysc(ix,jx)
     &         ,daz(ix,jx)

      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)
      dimension ss(ix,jx)

      dimension x(ix),xm(ix)
      dimension y(jx),ym(jx)
      dimension siny(jx),sinym(jx)

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
         drxsc(i,j) = 0.0
         drysc(i,j) = 0.0
         dbxsc(i,j) = 0.0
         dbysc(i,j) = 0.0
         daz(i,j) = 0.0
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
         pr(i,j) = ro(i,j)*cs2
      enddo
      enddo
      do j=1,jx
      do i=1,ix
         ez(i,j) = -vx(i,j)*by(i,j)+vy(i,j)*bx(i,j)
      enddo
      enddo
      do j=1,jx
      do i=1,ix
         rosc(i,j) = ro(i,j)*x(i)**2*siny(j)
         rxsc(i,j) = rx(i,j)*x(i)**2*siny(j)
         rysc(i,j) = ry(i,j)*x(i)*siny(j)
         bxsc(i,j) = bx(i,j)*x(i)**2*siny(j)
         bysc(i,j) = by(i,j)*x(i)*siny(j)
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

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)
     &            +pi8i*(by(i,j)**2-bx(i,j)**2)
         fy(i,j)= ro(i,j)*vx(i,j)*vy(i,j)-pi4i*bx(i,j)*by(i,j)
         fx(i,j)= fx(i,j)*x(i)**2*siny(j)
         fy(i,j)= fy(i,j)*x(i)*siny(j)
         ss(i,j)= (ro(i,j)*vy(i,j)**2+2*pr(i,j)
     &            +pi4i*(bx(i,j)**2))*x(i)*siny(j)
      enddo
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rxsch,drxsc,dt,ss,ix,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)-pi4i*by(i,j)*bx(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)**2+pr(i,j)
     &            +pi8i*(bx(i,j)**2-by(i,j)**2)
         fx(i,j)= fx(i,j)*x(i)*siny(j)
         fy(i,j)= fy(i,j)*siny(j)
         ss(i,j)= (ro(i,j)*vx(i,j)*vy(i,j)
     &            -pi4i*bx(i,j)*by(i,j))*(-2)*siny(j)
      enddo
      enddo
      call mlwhalf(rysc,rysch,drysc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rysch,drysc,dt,ss,ix,jx)

c---  x-magnetic ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= 0.
         fy(i,j)= ez(i,j)
         fy(i,j)= fy(i,j)*x(i)*siny(j)
      enddo
      enddo
      call mlwhalf(bxsc,bxsch,dbxsc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  y-magnetic ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= -ez(i,j)
         fy(i,j)= 0.
         fx(i,j)= fx(i,j)*x(i)*siny(j)
      enddo
      enddo
      call mlwhalf(bysc,bysch,dbysc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  z-magnetic potential ---
      do j=1,jx
      do i=1,ix
         ss(i,j)= -ez(i,j)
      enddo
      enddo
      call mlwsrch(azh ,daz ,dt,ss,ix,jx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         roh(i,j)   = rosch(i,j)/xm(i)**2/sinym(j)
         vxh(i,j)   = rxsch(i,j)/roh(i,j)/xm(i)**2/sinym(j)
         vyh(i,j)   = rysch(i,j)/roh(i,j)/xm(i)/sinym(j)
         bxh(i,j)   = bxsch(i,j)/xm(i)**2/sinym(j)
         byh(i,j)   = bysch(i,j)/xm(i)/sinym(j)
         prh(i,j)   = roh(i,j)*cs2
      enddo
      enddo
      do j=1,jx
      do i=1,ix
        ezh(i,j)=-vxh(i,j)*byh(i,j)+vyh(i,j)*bxh(i,j)
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

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+prh(i,j)
     &            +pi8i*(byh(i,j)**2-bxh(i,j)**2)
         fy(i,j)= roh(i,j)*vxh(i,j)*vyh(i,j)
     &            -pi4i*bxh(i,j)*byh(i,j)
         fx(i,j)= fx(i,j)*xm(i)**2*sinym(j)
         fy(i,j)= fy(i,j)*xm(i)*sinym(j)
         ss(i,j)= (roh(i,j)*vyh(i,j)**2+2*prh(i,j)
     &            +pi4i*(bxh(i,j)**2))*xm(i)*sinym(j)
      enddo
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vyh(i,j)*vxh(i,j)-pi4i*byh(i,j)*bxh(i,j)
         fy(i,j)= roh(i,j)*vyh(i,j)**2+prh(i,j)
     &            +pi8i*(bxh(i,j)**2-byh(i,j)**2)
         fx(i,j)= fx(i,j)*sinym(j)
         fy(i,j)= fy(i,j)*sinym(j)
         ss(i,j)= (roh(i,j)*vxh(i,j)*vyh(i,j)
     &            -pi4i*bxh(i,j)*byh(i,j))*(-2)*sinym(j)
      enddo
      enddo
      call mlwfull(drysc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drysc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  x-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= 0.
         fy(i,j)= ezh(i,j)
         fy(i,j)= fy(i,j)*xm(i)*sinym(j)
      enddo
      enddo
      call mlwfull(dbxsc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  y-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= -ezh(i,j)
         fy(i,j)= 0.
         fx(i,j)= fx(i,j)*xm(i)*sinym(j)
      enddo
      enddo
      call mlwfull(dbysc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  z-magnetic potential ---
      do j=1,jx
      do i=1,ix
         ss(i,j)= -ezh(i,j)
      enddo
      enddo
      call mlwsrcf(daz,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

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
      call mlwartv(rxsc,drxsc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rysc,drysc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(bxsc,dbxsc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(bysc,dbysc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         rosc(i,j) = rosc(i,j) +drosc(i,j)
         rxsc(i,j) = rxsc(i,j) +drxsc(i,j)
         rysc(i,j) = rysc(i,j) +drysc(i,j)
         bxsc(i,j) = bxsc(i,j) +dbxsc(i,j)
         bysc(i,j) = bysc(i,j) +dbysc(i,j)
         az(i,j)   = az(i,j)   +daz(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=2,jx-1  
      do i=2,ix-1  
         ro(i,j) = rosc(i,j)/x(i)**2/siny(j)
         vx(i,j) = rxsc(i,j)/ro(i,j)/x(i)**2/siny(j)
         vy(i,j) = rysc(i,j)/ro(i,j)/x(i)/siny(j)
         bx(i,j) = bxsc(i,j)/x(i)**2/siny(j)
         by(i,j) = bysc(i,j)/x(i)/siny(j)
      enddo
      enddo

      return
      end
