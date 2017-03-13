c======================================================================|
      subroutine mlw_m_s(ro,pr,vx,vy,vz,bx,by,bz,dt,qav,gm
     &     ,x,xm,dx,dxm,ix,y,ym,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  mlw_m_s
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * Spherical coordinate
c        
c INPUTS & OUTPUTS
c    ro(ix,jx,kx): [double] density
c    pr(ix,jx,kx): [double] pressure
c    vx(ix,jx,kx): [double] velocity
c    vy(ix,jx,kx): [double] velocity
c    vz(ix,jx,kx): [double] velocity
c    bx(ix,jx,kx): [double] magnetic field
c    by(ix,jx,kx): [double] magnetic field
c    bz(ix,jx,kx): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c    
c    x(ix), xm(ix): [double] coordinate (r)
c    y(jx), ym(jx): [double] coordinate (theta)
c    z(kx), zm(kx): [double] coordinate (phi)
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dz(kx), dzm(kx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
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
      dimension ro(ix,jx,kx),ee(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension roh(ix,jx,kx),eeh(ix,jx,kx),prh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)
      dimension exh(ix,jx,kx),eyh(ix,jx,kx),ezh(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)

      dimension rosc(ix,jx,kx),eesc(ix,jx,kx)
     &         ,rxsc(ix,jx,kx),rysc(ix,jx,kx),rzsc(ix,jx,kx)
     &         ,bxsc(ix,jx,kx),bysc(ix,jx,kx),bzsc(ix,jx,kx)
      dimension rosch(ix,jx,kx),eesch(ix,jx,kx)
     &         ,rxsch(ix,jx,kx),rysch(ix,jx,kx),rzsch(ix,jx,kx)
     &         ,bxsch(ix,jx,kx),bysch(ix,jx,kx),bzsch(ix,jx,kx)
      dimension drosc(ix,jx,kx),deesc(ix,jx,kx)
     &         ,drxsc(ix,jx,kx),drysc(ix,jx,kx),drzsc(ix,jx,kx)
     &         ,dbxsc(ix,jx,kx),dbysc(ix,jx,kx),dbzsc(ix,jx,kx)
      dimension fx(ix,jx,kx),qx(ix,jx,kx)
      dimension fy(ix,jx,kx),qy(ix,jx,kx)
      dimension fz(ix,jx,kx),qz(ix,jx,kx)
      dimension ss(ix,jx,kx)
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
         deesc(i,j,k)= 0.0
         drxsc(i,j,k)= 0.0
         drysc(i,j,k)= 0.0
         drzsc(i,j,k)= 0.0
         dbxsc(i,j,k)= 0.0
         dbysc(i,j,k)= 0.0
         dbzsc(i,j,k)= 0.0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         vv=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         bb=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         ee(i,j,k) = pr(i,j,k)/(gm-1)+0.5*ro(i,j,k)*vv+pi8i*bb
         rx(i,j,k) = ro(i,j,k)*vx(i,j,k)
         ry(i,j,k) = ro(i,j,k)*vy(i,j,k)
         rz(i,j,k) = ro(i,j,k)*vz(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ex(i,j,k) = -vy(i,j,k)*bz(i,j,k)+vz(i,j,k)*by(i,j,k)
         ey(i,j,k) = -vz(i,j,k)*bx(i,j,k)+vx(i,j,k)*bz(i,j,k)
         ez(i,j,k) = -vx(i,j,k)*by(i,j,k)+vy(i,j,k)*bx(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
         rosc(i,j,k) = ro(i,j,k)*x(i)**2*siny(j)
         eesc(i,j,k) = ee(i,j,k)*x(i)**2*siny(j)
         rxsc(i,j,k) = rx(i,j,k)*x(i)**2*siny(j)
         rysc(i,j,k) = ry(i,j,k)*x(i)*siny(j)
         rzsc(i,j,k) = rz(i,j,k)*x(i)
         bxsc(i,j,k) = bx(i,j,k)*x(i)**2*siny(j)
         bysc(i,j,k) = by(i,j,k)*x(i)*siny(j)
         bzsc(i,j,k) = bz(i,j,k)*x(i)
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

c---  energy ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         vv=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         ep= pr(i,j,k)*gm/(gm-1.)+0.5*ro(i,j,k)*vv
         px= (bz(i,j,k)*ey(i,j,k)-by(i,j,k)*ez(i,j,k))*pi4i
         py= (bx(i,j,k)*ez(i,j,k)-bz(i,j,k)*ex(i,j,k))*pi4i
         pz= (by(i,j,k)*ex(i,j,k)-bx(i,j,k)*ey(i,j,k))*pi4i
         fx(i,j,k)= ep*vx(i,j,k)+px
         fy(i,j,k)= ep*vy(i,j,k)+py
         fz(i,j,k)= ep*vz(i,j,k)+pz
         fx(i,j,k)= fx(i,j,k)*x(i)**2*siny(j)
         fy(i,j,k)= fy(i,j,k)*x(i)*siny(j)
         fz(i,j,k)= fz(i,j,k)*x(i)
      enddo
      enddo
      enddo
      call mlwhalf(eesc,eesch,deesc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
     &          +pi8i*(by(i,j,k)**2+bz(i,j,k)**2-bx(i,j,k)**2)
         fy(i,j,k)= ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
     &          -pi4i*bx(i,j,k)*by(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
     &          -pi4i*bx(i,j,k)*bz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)**2*siny(j)
         fy(i,j,k)= fy(i,j,k)*x(i)*siny(j)
         fz(i,j,k)= fz(i,j,k)*x(i)
         ss(i,j,k)= 
     &            (ro(i,j,k)*(vy(i,j,k)**2+vz(i,j,k)**2)+2*pr(i,j,k)
     &            +pi4i*bx(i,j,k)**2)*x(i)*siny(j)
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
     &          -pi4i*by(i,j,k)*bx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)**2+pr(i,j,k)
     &          +pi8i*(bx(i,j,k)**2+bz(i,j,k)**2-by(i,j,k)**2)
         fz(i,j,k)= ro(i,j,k)*vy(i,j,k)*vz(i,j,k)
     &          -pi4i*by(i,j,k)*bz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*x(i)*siny(j)
         fy(i,j,k)= fy(i,j,k)*siny(j)
         ss(i,j,k)= 
     &            (ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
     &            -pi4i*bx(i,j,k)*by(i,j,k))*(-2)*siny(j)
     &          + (ro(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
     &            +pi8i*(bx(i,j,k)**2+by(i,j,k)**2-bz(i,j,k)**2))
     &            *cos(y(j))
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
     &          -pi4i*bz(i,j,k)*bx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vz(i,j,k)*vy(i,j,k)
     &          -pi4i*bz(i,j,k)*by(i,j,k)
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
     &          +pi8i*(bx(i,j,k)**2+by(i,j,k)**2-bz(i,j,k)**2)
         fx(i,j,k)= fx(i,j,k)*x(i)
         fz(i,j,k)= fz(i,j,k)/siny(j)
         ss(i,j,k)=
     &            (ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
     &            -pi4i*bx(i,j,k)*bz(i,j,k))*(-2)
     &          + (ro(i,j,k)*vy(i,j,k)*vz(i,j,k)
     &            -pi4i*by(i,j,k)*bz(i,j,k))*(-2)/tan(y(j))
      enddo
      enddo
      enddo
      call mlwhalf(rzsc,rzsch,drzsc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rzsch,drzsc,dt,ss,ix,jx,kx)

c---  x-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  0.
         fy(i,j,k)=  ez(i,j,k)*x(i)*siny(j)
         fz(i,j,k)= -ey(i,j,k)*x(i)
      enddo
      enddo
      enddo
      call mlwhalf(bxsc,bxsch,dbxsc,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  y-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= -ez(i,j,k)*x(i)*siny(j)
         fy(i,j,k)=  0.
         fz(i,j,k)=  ex(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(bysc,bysch,dbysc,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  z-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  ey(i,j,k)*x(i)
         fy(i,j,k)= -ex(i,j,k)
         fz(i,j,k)=  0.
      enddo
      enddo
      enddo
      call mlwhalf(bzsc,bzsch,dbzsc,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         roh(i,j,k) = rosch(i,j,k)/xm(i)**2/sinym(j)
         eeh(i,j,k) = eesch(i,j,k)/xm(i)**2/sinym(j)
         vxh(i,j,k) = rxsch(i,j,k)/roh(i,j,k)/xm(i)**2/sinym(j)
         vyh(i,j,k) = rysch(i,j,k)/roh(i,j,k)/xm(i)/sinym(j)
         vzh(i,j,k) = rzsch(i,j,k)/roh(i,j,k)/xm(i)
         bxh(i,j,k) = bxsch(i,j,k)/xm(i)**2/sinym(j)
         byh(i,j,k) = bysch(i,j,k)/xm(i)/sinym(j)
         bzh(i,j,k) = bzsch(i,j,k)/xm(i)
         vv=vxh(i,j,k)**2+vyh(i,j,k)**2+vzh(i,j,k)**2
         bb=bxh(i,j,k)**2+byh(i,j,k)**2+bzh(i,j,k)**2
         prh(i,j,k)   = (gm-1)*(eeh(i,j,k)-0.5*roh(i,j,k)*vv-pi8i*bb)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
        exh(i,j,k)=-vyh(i,j,k)*bzh(i,j,k)+vzh(i,j,k)*byh(i,j,k)
        eyh(i,j,k)=-vzh(i,j,k)*bxh(i,j,k)+vxh(i,j,k)*bzh(i,j,k)
        ezh(i,j,k)=-vxh(i,j,k)*byh(i,j,k)+vyh(i,j,k)*bxh(i,j,k)
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

c---  energy     ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         vv=vxh(i,j,k)**2+vyh(i,j,k)**2+vzh(i,j,k)**2
         ep= prh(i,j,k)*gm/(gm-1.)+0.5*roh(i,j,k)*vv
         px= (bzh(i,j,k)*eyh(i,j,k)-byh(i,j,k)*ezh(i,j,k))*pi4i
         py= (bxh(i,j,k)*ezh(i,j,k)-bzh(i,j,k)*exh(i,j,k))*pi4i
         pz= (byh(i,j,k)*exh(i,j,k)-bxh(i,j,k)*eyh(i,j,k))*pi4i
         fx(i,j,k)= ep*vxh(i,j,k)+px
         fy(i,j,k)= ep*vyh(i,j,k)+py
         fz(i,j,k)= ep*vzh(i,j,k)+pz
         fx(i,j,k)= fx(i,j,k)*xm(i)**2*sinym(j)
         fy(i,j,k)= fy(i,j,k)*xm(i)*sinym(j)
         fz(i,j,k)= fz(i,j,k)*xm(i)
      enddo
      enddo
      enddo
      call mlwfull(deesc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= roh(i,j,k)*vxh(i,j,k)**2+prh(i,j,k)
     &          +pi8i*(byh(i,j,k)**2+bzh(i,j,k)**2-bxh(i,j,k)**2)
         fy(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
     &          -pi4i*bxh(i,j,k)*byh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
     &          -pi4i*bxh(i,j,k)*bzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)**2*sinym(j)
         fy(i,j,k)= fy(i,j,k)*xm(i)*sinym(j)
         fz(i,j,k)= fz(i,j,k)*xm(i)
         ss(i,j,k)=
     &            (roh(i,j,k)*(vyh(i,j,k)**2+vzh(i,j,k)**2)
     &            +2*prh(i,j,k)
     &            +pi4i*bxh(i,j,k)**2)*xm(i)*sinym(j)
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
     &          -pi4i*byh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vyh(i,j,k)**2+prh(i,j,k)
     &          +pi8i*(bxh(i,j,k)**2+bzh(i,j,k)**2-byh(i,j,k)**2)
         fz(i,j,k)= roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
     &          -pi4i*byh(i,j,k)*bzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*xm(i)*sinym(j)
         fy(i,j,k)= fy(i,j,k)*sinym(j)
         ss(i,j,k)=
     &            (roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
     &            -pi4i*bxh(i,j,k)*byh(i,j,k))*(-2)*sinym(j)
     &          + (roh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
     &            +pi8i*(bxh(i,j,k)**2+byh(i,j,k)**2-bzh(i,j,k)**2))
     &            *cos(ym(j))
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
     &          -pi4i*bzh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vzh(i,j,k)*vyh(i,j,k)
     &          -pi4i*bzh(i,j,k)*byh(i,j,k)
         fz(i,j,k)= roh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
     &          +pi8i*(bxh(i,j,k)**2+byh(i,j,k)**2-bzh(i,j,k)**2)
         fx(i,j,k)= fx(i,j,k)*xm(i)
         fz(i,j,k)= fz(i,j,k)/sinym(j)
         ss(i,j,k)=
     &            (roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
     &            -pi4i*bxh(i,j,k)*bzh(i,j,k))*(-2)
     &          + (roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
     &            -pi4i*byh(i,j,k)*bzh(i,j,k))*(-2)/tan(ym(j))
      enddo
      enddo
      enddo
      call mlwfull(drzsc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drzsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  x-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  0.
         fy(i,j,k)=  ezh(i,j,k)*xm(i)*sinym(j)
         fz(i,j,k)= -eyh(i,j,k)*xm(i)
      enddo
      enddo
      enddo
      call mlwfull(dbxsc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  y-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= -ezh(i,j,k)*xm(i)*sinym(j)
         fy(i,j,k)=  0.
         fz(i,j,k)=  exh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dbysc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  z-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  eyh(i,j,k)*xm(i)
         fy(i,j,k)= -exh(i,j,k)
         fz(i,j,k)=  0.
      enddo
      enddo
      enddo
      call mlwfull(dbzsc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

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
      call mlwartv(eesc,deesc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rxsc,drxsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rysc,drysc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rzsc,drzsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(bxsc,dbxsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(bysc,dbysc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(bzsc,dbzsc,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         rosc(i,j,k)= rosc(i,j,k) +drosc(i,j,k)
         eesc(i,j,k)= eesc(i,j,k) +deesc(i,j,k)
         rxsc(i,j,k)= rxsc(i,j,k) +drxsc(i,j,k)
         rysc(i,j,k)= rysc(i,j,k) +drysc(i,j,k)
         rzsc(i,j,k)= rzsc(i,j,k) +drzsc(i,j,k)
         bxsc(i,j,k)= bxsc(i,j,k) +dbxsc(i,j,k)
         bysc(i,j,k)= bysc(i,j,k) +dbysc(i,j,k)
         bzsc(i,j,k)= bzsc(i,j,k) +dbzsc(i,j,k)
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
         ee(i,j,k) = eesc(i,j,k)/x(i)**2/siny(j)
         vx(i,j,k) = rxsc(i,j,k)/ro(i,j,k)/x(i)**2/siny(j)
         vy(i,j,k) = rysc(i,j,k)/ro(i,j,k)/x(i)/siny(j)
         vz(i,j,k) = rzsc(i,j,k)/ro(i,j,k)/x(i)
         bx(i,j,k) = bxsc(i,j,k)/x(i)**2/siny(j)
         by(i,j,k) = bysc(i,j,k)/x(i)/siny(j)
         bz(i,j,k) = bzsc(i,j,k)/x(i)
         vv=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         bb=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         pr(i,j,k) = (gm-1)*(ee(i,j,k) - 0.5*ro(i,j,k)*vv - pi8i*bb)
      enddo
      enddo
      enddo

      return
      end
