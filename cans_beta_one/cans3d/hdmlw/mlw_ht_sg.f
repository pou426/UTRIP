c======================================================================|
      subroutine mlw_ht_sg(ro,vx,vy,vz,dt,qav,cs2
     &     ,gx,gxm,gy,gym,gz,gzm
     &     ,x,xm,dx,dxm,ix,y,ym,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  mlw_ht_sg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal hydrodynamics
c        * Spherical coordinate
c        * gravity
c        
c INPUTS & OUTPUTS
c    ro(ix,jx,kx): [double] density
c    vx(ix,jx,kx): [double] velocity
c    vy(ix,jx,kx): [double] velocity
c    vz(ix,jx,kx): [double] velocity
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c    
c    gx(ix,jx,kx), gxm(ix,jx,kx) : [double] gravity
c    gy(ix,jx,kx), gym(ix,jx,kx) : [double] gravity
c    gz(ix,jx,kx), gzm(ix,jx,kx) : [double] gravity
c    x(ix), xm(ix): [double] coordinate (r)
c    y(jx), ym(jx): [double] coordinate (theta)
c    z(kx), zm(kx): [double] coordinate (phi)
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dz(kx), dzm(kx): [double] grid spacing
c    dt: [double] delta time
c    cs2: [double] square of sound speed
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension roh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)

      dimension rosc(ix,jx,kx)
     &         ,rxsc(ix,jx,kx),rysc(ix,jx,kx),rzsc(ix,jx,kx)
      dimension rosch(ix,jx,kx)
     &         ,rxsch(ix,jx,kx),rysch(ix,jx,kx),rzsch(ix,jx,kx)
      dimension drosc(ix,jx,kx)
     &         ,drxsc(ix,jx,kx),drysc(ix,jx,kx),drzsc(ix,jx,kx)
      dimension fx(ix,jx,kx),qx(ix,jx,kx)
      dimension fy(ix,jx,kx),qy(ix,jx,kx)
      dimension fz(ix,jx,kx),qz(ix,jx,kx)
      dimension ss(ix,jx,kx)
      dimension x(ix),xm(ix)
      dimension y(jx),ym(jx)
      dimension siny(jx),sinym(jx)
      dimension gx(ix,jx,kx),gxm(ix,jx,kx)
      dimension gy(ix,jx,kx),gym(ix,jx,kx)
      dimension gz(ix,jx,kx),gzm(ix,jx,kx)
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

      do k=1,kx
         dzi(k)  = 1.0/dz(k)
         dzim(k) = 1.0/dzm(k)
      enddo
      do k=2,kx-1
         uz1(k)  = 0.5*dzm(k-1)/dz(k)
         uz0(k)  = 0.5*dzm(k)/dz(k)
      enddo

      do j=1,jx
         siny(j) = sin(y(j))
         sinym(j) = sin(ym(j))
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         drosc(i,j,k)= 0.0
         drxsc(i,j,k)= 0.0
         drysc(i,j,k)= 0.0
         drzsc(i,j,k)= 0.0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         rx(i,j,k) = ro(i,j,k)*vx(i,j,k)
         ry(i,j,k) = ro(i,j,k)*vy(i,j,k)
         rz(i,j,k) = ro(i,j,k)*vz(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
         rosc(i,j,k) = ro(i,j,k)*x(i)**2*siny(j)
         rxsc(i,j,k) = rx(i,j,k)*x(i)**2*siny(j)
         rysc(i,j,k) = ry(i,j,k)*x(i)*siny(j)
         rzsc(i,j,k) = rz(i,j,k)*x(i)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)**2*siny(j)
         fy(i,j,k)= fy(i,j,k)*x(i)*siny(j)
         fz(i,j,k)= fz(i,j,k)*x(i)
      enddo
      enddo
      enddo
      call mlwhalf(rosc,rosch,drosc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)**2+cs2*ro(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)**2*siny(j)
         fy(i,j,k)= fy(i,j,k)*x(i)*siny(j)
         fz(i,j,k)= fz(i,j,k)*x(i)
         ss(i,j,k)= 
     &        (ro(i,j,k)*(vy(i,j,k)**2+vz(i,j,k)**2)+2*cs2*ro(i,j,k)
     &                              )*x(i)*siny(j)
     &          +  ro(i,j,k)*gx(i,j,k)*x(i)**2*siny(j)
      enddo
      enddo
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rxsch,drxsc,dt,ss,ix,jx,kx)

c---  y-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vy(i,j,k)*vx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)**2+cs2*ro(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vy(i,j,k)*vz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)*siny(j)
         fy(i,j,k)= fy(i,j,k)*siny(j)
         ss(i,j,k)= 
     &            (ro(i,j,k)*vx(i,j,k)*vy(i,j,k))*(-2)*siny(j)
     &          + (ro(i,j,k)*vz(i,j,k)**2+cs2*ro(i,j,k))*cos(y(j))
     &          +  ro(i,j,k)*gy(i,j,k)*x(i)*siny(j)
      enddo
      enddo
      enddo
      call mlwhalf(rysc,rysch,drysc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rysch,drysc,dt,ss,ix,jx,kx)

c---  z-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vz(i,j,k)*vx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vz(i,j,k)*vy(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)**2+cs2*ro(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)
         fz(i,j,k)= fz(i,j,k)/siny(j)
         ss(i,j,k)=
     &            (ro(i,j,k)*vx(i,j,k)*vz(i,j,k))*(-2)
     &          + (ro(i,j,k)*vy(i,j,k)*vz(i,j,k))*(-2)/tan(y(j))
     &          +  ro(i,j,k)*gz(i,j,k)*x(i)
      enddo
      enddo
      enddo
      call mlwhalf(rzsc,rzsch,drzsc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rzsch,drzsc,dt,ss,ix,jx,kx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         roh(i,j,k) = rosch(i,j,k)/xm(i)**2/sinym(j)
         vxh(i,j,k) = rxsch(i,j,k)/roh(i,j,k)/xm(i)**2/sinym(j)
         vyh(i,j,k) = rysch(i,j,k)/roh(i,j,k)/xm(i)/sinym(j)
         vzh(i,j,k) = rzsch(i,j,k)/roh(i,j,k)/xm(i)
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= roh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vyh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)**2*sinym(j)
         fy(i,j,k)= fy(i,j,k)*xm(i)*sinym(j)
         fz(i,j,k)= fz(i,j,k)*xm(i)
      enddo
      enddo
      enddo
      call mlwfull(drosc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= roh(i,j,k)*vxh(i,j,k)**2+cs2*roh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)**2*sinym(j)
         fy(i,j,k)= fy(i,j,k)*xm(i)*sinym(j)
         fz(i,j,k)= fz(i,j,k)*xm(i)
         ss(i,j,k)=
     &            (roh(i,j,k)*(vyh(i,j,k)**2+vzh(i,j,k)**2)
     &            +2*cs2*roh(i,j,k)
     &                               )*xm(i)*sinym(j)
     &          +  roh(i,j,k)*gxm(i,j,k)*xm(i)**2*sinym(j)
      enddo
      enddo
      enddo
      call mlwfull(drxsc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  y-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= roh(i,j,k)*vyh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vyh(i,j,k)**2+cs2*roh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)*sinym(j)
         fy(i,j,k)= fy(i,j,k)*sinym(j)
         ss(i,j,k)=
     &            (roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k))*(-2)*sinym(j)
     &       + (roh(i,j,k)*vzh(i,j,k)**2+cs2*roh(i,j,k))*cos(ym(j))
     &          +  roh(i,j,k)*gym(i,j,k)*xm(i)*sinym(j)
      enddo
      enddo
      enddo
      call mlwfull(drysc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drysc,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  z-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= roh(i,j,k)*vzh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vzh(i,j,k)*vyh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vzh(i,j,k)**2+cs2*roh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)
         fz(i,j,k)= fz(i,j,k)/sinym(j)
         ss(i,j,k)=
     &            (roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k))*(-2)
     &          + (roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k))*(-2)/tan(ym(j))
     &          +  roh(i,j,k)*gzm(i,j,k)*xm(i)
      enddo
      enddo
      enddo
      call mlwfull(drzsc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drzsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=1.0
      zero=0.0d0
      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        qx(i,j,k)=qav*dxm(i)*max(zero,abs(vx(i+1,j,k)-vx(i,j,k))-1.0e-4)
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        qy(i,j,k)=qav*dym(j)*max(zero,abs(vy(i,j+1,k)-vy(i,j,k))-1.0e-4)
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        qz(i,j,k)=qav*dzm(k)*max(zero,abs(vz(i,j,k+1)-vz(i,j,k))-1.0e-4)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(rosc,drosc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rxsc,drxsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rysc,drysc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rzsc,drzsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         rosc(i,j,k)= rosc(i,j,k) +drosc(i,j,k)
         rxsc(i,j,k)= rxsc(i,j,k) +drxsc(i,j,k)
         rysc(i,j,k)= rysc(i,j,k) +drysc(i,j,k)
         rzsc(i,j,k)= rzsc(i,j,k) +drzsc(i,j,k)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do k=2,kx-1  
      do j=2,jx-1  
      do i=2,ix-1  
         ro(i,j,k) = rosc(i,j,k)/x(i)**2/siny(j)
         vx(i,j,k) = rxsc(i,j,k)/ro(i,j,k)/x(i)**2/siny(j)
         vy(i,j,k) = rysc(i,j,k)/ro(i,j,k)/x(i)/siny(j)
         vz(i,j,k) = rzsc(i,j,k)/ro(i,j,k)/x(i)
      enddo
      enddo
      enddo

      return
      end
