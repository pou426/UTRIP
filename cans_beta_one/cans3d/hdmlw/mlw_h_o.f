c======================================================================|
      subroutine mlw_h_o(ro,pr,vx,vy,vz,dt,qav,gm
     &                  ,hh,hhm,gi,gim,gd,gdm
     &                  ,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  mlw_h_o
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * Hydrodynamics
c        * orthogonal general coordinate
c        
c INPUTS & OUTPUTS
c    ro(ix,jx,kx): [double] density
c    pr(ix,jx,kx): [double] pressure
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
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dz(kx), dzm(kx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix,jx,kx: [integer] dimension size
c    
c HISTORY
c    written 2004-3-15 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension roh(ix,jx,kx),prh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)

      dimension rosc(ix,jx,kx),eesc(ix,jx,kx)
     &         ,rxsc(ix,jx,kx),rysc(ix,jx,kx),rzsc(ix,jx,kx)
      dimension rosch(ix,jx,kx),eesch(ix,jx,kx)
     &         ,rxsch(ix,jx,kx),rysch(ix,jx,kx),rzsch(ix,jx,kx)
      dimension drosc(ix,jx,kx),deesc(ix,jx,kx)
     &         ,drxsc(ix,jx,kx),drysc(ix,jx,kx),drzsc(ix,jx,kx)
      dimension fx(ix,jx,kx),qx(ix,jx,kx)
      dimension fy(ix,jx,kx),qy(ix,jx,kx)
      dimension fz(ix,jx,kx),qz(ix,jx,kx)
      dimension ss(ix,jx,kx)
      dimension hh(ix,jx,kx,3),gi(ix,jx,kx,3),gd(ix,jx,kx,3)
      dimension hhm(ix,jx,kx,3),gim(ix,jx,kx,3),gdm(ix,jx,kx,3)
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
c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         drosc(i,j,k)= 0.d0
         deesc(i,j,k)= 0.d0
         drxsc(i,j,k)= 0.d0
         drysc(i,j,k)= 0.d0
         drzsc(i,j,k)= 0.d0
         rosc(i,j,k) = 0.d0
         eesc(i,j,k) = 0.d0
         rxsc(i,j,k) = 0.d0
         rysc(i,j,k) = 0.d0
         rzsc(i,j,k) = 0.d0
         rosch(i,j,k)= 0.d0
         eesch(i,j,k)= 0.d0
         rxsch(i,j,k)= 0.d0
         rysch(i,j,k)= 0.d0
         rzsch(i,j,k)= 0.d0
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
         ee        = pr(i,j,k)/(gm-1)+0.5*ro(i,j,k)*vv
         rx        = ro(i,j,k)*vx(i,j,k)
         ry        = ro(i,j,k)*vy(i,j,k)
         rz        = ro(i,j,k)*vz(i,j,k)

         sc=hh(i,j,k,1)*hh(i,j,k,2)*hh(i,j,k,3)
         rosc(i,j,k) = ro(i,j,k)*sc
         eesc(i,j,k) = ee*sc
         rxsc(i,j,k) = rx*sc
         rysc(i,j,k) = ry*sc
         rzsc(i,j,k) = rz*sc
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
         fx(i,j,k)= fx(i,j,k)*hh(i,j,k,2)*hh(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hh(i,j,k,3)*hh(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)
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
         fx(i,j,k)= ep*vx(i,j,k)
         fy(i,j,k)= ep*vy(i,j,k)
         fz(i,j,k)= ep*vz(i,j,k)
         fx(i,j,k)= fx(i,j,k)*hh(i,j,k,2)*hh(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hh(i,j,k,3)*hh(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)
      enddo
      enddo
      enddo
      call mlwhalf(eesc,eesch,deesc,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         t11= ro(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
         t22= ro(i,j,k)*vy(i,j,k)**2+pr(i,j,k)
         t33= ro(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
         t12= ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
         t13= ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
         t21= t12
         t31= t13
         fx(i,j,k)= t11
         fy(i,j,k)= t21
         fz(i,j,k)= t31
         fx(i,j,k)= fx(i,j,k)*hh(i,j,k,2)*hh(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hh(i,j,k,3)*hh(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)
         ss(i,j,k)=-gd(i,j,k,3)*t22-gi(i,j,k,2)*t33
     &             +gi(i,j,k,3)*t12+gd(i,j,k,2)*t13
         ss(i,j,k)= ss(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)*hh(i,j,k,3)
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
         t11= ro(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
         t22= ro(i,j,k)*vy(i,j,k)**2+pr(i,j,k)
         t33= ro(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
         t23= ro(i,j,k)*vy(i,j,k)*vz(i,j,k)
         t21= ro(i,j,k)*vy(i,j,k)*vx(i,j,k)
         t32= t23
         t12= t21
         fx(i,j,k)= t12
         fy(i,j,k)= t22
         fz(i,j,k)= t32
         fx(i,j,k)= fx(i,j,k)*hh(i,j,k,2)*hh(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hh(i,j,k,3)*hh(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)
         ss(i,j,k)=-gd(i,j,k,1)*t33-gi(i,j,k,3)*t11
     &             +gi(i,j,k,1)*t23+gd(i,j,k,3)*t21
         ss(i,j,k)= ss(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)*hh(i,j,k,3)
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
         t11= ro(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
         t22= ro(i,j,k)*vy(i,j,k)**2+pr(i,j,k)
         t33= ro(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
         t32= ro(i,j,k)*vz(i,j,k)*vy(i,j,k)
         t31= ro(i,j,k)*vz(i,j,k)*vx(i,j,k)
         t23= t32
         t13= t31
         fx(i,j,k)= t13
         fy(i,j,k)= t23
         fz(i,j,k)= t33
         fx(i,j,k)= fx(i,j,k)*hh(i,j,k,2)*hh(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hh(i,j,k,3)*hh(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)
         ss(i,j,k)=-gd(i,j,k,2)*t33-gi(i,j,k,1)*t11
     &             +gi(i,j,k,2)*t23+gd(i,j,k,1)*t21
         ss(i,j,k)= ss(i,j,k)*hh(i,j,k,1)*hh(i,j,k,2)*hh(i,j,k,3)
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
         scm=hhm(i,j,k,1)*hhm(i,j,k,2)*hhm(i,j,k,3)
         roh(i,j,k) = rosch(i,j,k)/scm
         eeh        = eesch(i,j,k)/scm
         rxh        = rxsch(i,j,k)/scm
         ryh        = rysch(i,j,k)/scm
         rzh        = rzsch(i,j,k)/scm
         vxh(i,j,k) = rxh/roh(i,j,k)
         vyh(i,j,k) = ryh/roh(i,j,k)
         vzh(i,j,k) = rzh/roh(i,j,k)
         vv=vxh(i,j,k)**2+vyh(i,j,k)**2+vzh(i,j,k)**2
         prh(i,j,k)   = (gm-1)*(eeh-0.5*roh(i,j,k)*vv)
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
         fx(i,j,k)= fx(i,j,k)*hhm(i,j,k,2)*hhm(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hhm(i,j,k,3)*hhm(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)
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
         fx(i,j,k)= ep*vxh(i,j,k)
         fy(i,j,k)= ep*vyh(i,j,k)
         fz(i,j,k)= ep*vzh(i,j,k)
         fx(i,j,k)= fx(i,j,k)*hhm(i,j,k,2)*hhm(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hhm(i,j,k,3)*hhm(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)
      enddo
      enddo
      enddo
      call mlwfull(deesc,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         t11= roh(i,j,k)*vxh(i,j,k)**2+prh(i,j,k)
         t22= roh(i,j,k)*vyh(i,j,k)**2+prh(i,j,k)
         t33= roh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
         t12= roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
         t13= roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
         t21= t12
         t31= t13
         fx(i,j,k)= t11
         fy(i,j,k)= t21
         fz(i,j,k)= t31
         fx(i,j,k)= fx(i,j,k)*hhm(i,j,k,2)*hhm(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hhm(i,j,k,3)*hhm(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)
         ss(i,j,k)=-gdm(i,j,k,3)*t22-gim(i,j,k,2)*t33
     &             +gim(i,j,k,3)*t12+gdm(i,j,k,2)*t13
         ss(i,j,k)= ss(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)*hhm(i,j,k,3)
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
         t11= roh(i,j,k)*vxh(i,j,k)**2+prh(i,j,k)
         t22= roh(i,j,k)*vyh(i,j,k)**2+prh(i,j,k)
         t33= roh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
         t23= roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
         t21= roh(i,j,k)*vyh(i,j,k)*vxh(i,j,k)
         t32= t23
         t12= t21
         fx(i,j,k)= t12
         fy(i,j,k)= t22
         fz(i,j,k)= t32
         fx(i,j,k)= fx(i,j,k)*hhm(i,j,k,2)*hhm(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hhm(i,j,k,3)*hhm(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)
         ss(i,j,k)=-gdm(i,j,k,1)*t33-gim(i,j,k,3)*t11
     &             +gim(i,j,k,1)*t23+gdm(i,j,k,3)*t21
         ss(i,j,k)= ss(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)*hhm(i,j,k,3)
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
         t11= roh(i,j,k)*vxh(i,j,k)**2+prh(i,j,k)
         t22= roh(i,j,k)*vyh(i,j,k)**2+prh(i,j,k)
         t33= roh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
         t32= roh(i,j,k)*vzh(i,j,k)*vyh(i,j,k)
         t31= roh(i,j,k)*vzh(i,j,k)*vxh(i,j,k)
         t23= t32
         t13= t31
         fx(i,j,k)= t13
         fy(i,j,k)= t23
         fz(i,j,k)= t33
         fx(i,j,k)= fx(i,j,k)*hhm(i,j,k,2)*hhm(i,j,k,3)
         fy(i,j,k)= fy(i,j,k)*hhm(i,j,k,3)*hhm(i,j,k,1)
         fz(i,j,k)= fz(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)
         ss(i,j,k)=-gdm(i,j,k,2)*t33-gim(i,j,k,1)*t11
     &             +gim(i,j,k,2)*t23+gdm(i,j,k,1)*t21
         ss(i,j,k)= ss(i,j,k)*hhm(i,j,k,1)*hhm(i,j,k,2)*hhm(i,j,k,3)
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
      call mlwartv(eesc,deesc,dt
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
         eesc(i,j,k)= eesc(i,j,k) +deesc(i,j,k)
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
         sc=hh(i,j,k,1)*hh(i,j,k,2)*hh(i,j,k,3)
         ro(i,j,k) = rosc(i,j,k)/sc
         ee        = eesc(i,j,k)/sc
         rx        = rxsc(i,j,k)/sc
         ry        = rysc(i,j,k)/sc
         rz        = rzsc(i,j,k)/sc
         vx(i,j,k) = rx/ro(i,j,k)
         vy(i,j,k) = ry/ro(i,j,k)
         vz(i,j,k) = rz/ro(i,j,k)
         vv=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         pr(i,j,k) = (gm-1)*(ee - 0.5*ro(i,j,k)*vv )
      enddo
      enddo
      enddo

      return
      end
