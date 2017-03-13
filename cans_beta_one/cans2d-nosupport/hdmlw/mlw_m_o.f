c======================================================================|
      subroutine mlw_m_o(ro,pr,vx,vy,bx,by,az,dt,qav,gm
     &                  ,hh,hhm,gi,gim,gd,gdm
     &                  ,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  mlw_m_o
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * orthogonal general coordinate
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
c    x(ix), xm(ix): [double] coordinate
c    y(jx), ym(jx): [double] coordinate
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2004-3-15 T. Yokoyama
c 
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
     &         ,bx(ix,jx),by(ix,jx)
      dimension ez(ix,jx),pm(ix,jx)
      dimension az(ix,jx)
      dimension prh(ix,jx),roh(ix,jx),vxh(ix,jx),vyh(ix,jx)
     &         ,bxh(ix,jx),byh(ix,jx),ezh(ix,jx),pmh(ix,jx)
     &         ,azh(ix,jx)

      dimension rosc(ix,jx),eesc(ix,jx),rxsc(ix,jx),rysc(ix,jx)
     &         ,bxsc(ix,jx),bysc(ix,jx)
      dimension rosch(ix,jx),eesch(ix,jx),rxsch(ix,jx),rysch(ix,jx)
     &         ,bxsch(ix,jx),bysch(ix,jx)
      dimension drosc(ix,jx),deesc(ix,jx),drxsc(ix,jx),drysc(ix,jx)
     &         ,dbxsc(ix,jx),dbysc(ix,jx)
     &         ,daz(ix,jx)

      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)
      dimension ss(ix,jx)
      dimension hh(ix,jx,3),gi(ix,jx,3),gd(ix,jx,3)
      dimension hhm(ix,jx,3),gim(ix,jx,3),gdm(ix,jx,3)
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
c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         drosc(i,j) = 0.d0
         deesc(i,j) = 0.d0
         drxsc(i,j) = 0.d0
         drysc(i,j) = 0.d0
         dbxsc(i,j) = 0.d0
         dbysc(i,j) = 0.d0
         daz(i,j)   = 0.d0
         rosc(i,j)  = 0.d0
         eesc(i,j)  = 0.d0
         rxsc(i,j)  = 0.d0
         rysc(i,j)  = 0.d0
         bxsc(i,j)  = 0.d0
         bysc(i,j)  = 0.d0
         rosch(i,j) = 0.d0
         eesch(i,j) = 0.d0
         rxsch(i,j) = 0.d0
         rysch(i,j) = 0.d0
         bxsch(i,j) = 0.d0
         bysch(i,j) = 0.d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vy(i,j)**2
         bb=bx(i,j)**2+by(i,j)**2
         pm(i,j) = pi8i*bb
         ee      = pr(i,j)/(gm-1)+0.5*ro(i,j)*vv+pm(i,j)
         rx      = ro(i,j)*vx(i,j)
         ry      = ro(i,j)*vy(i,j)
         ez(i,j) = -vx(i,j)*by(i,j)+vy(i,j)*bx(i,j)

         sc=hh(i,j,1)*hh(i,j,2)*hh(i,j,3)
         rosc(i,j) = ro(i,j)*sc
         eesc(i,j) = ee*sc
         rxsc(i,j) = rx*sc
         rysc(i,j) = ry*sc
         bxsc(i,j) = bx(i,j)*hh(i,j,2)*hh(i,j,3)
         bysc(i,j) = by(i,j)*hh(i,j,3)*hh(i,j,1)
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
         fx(i,j)= fx(i,j)*hh(i,j,2)*hh(i,j,3)
         fy(i,j)= fy(i,j)*hh(i,j,3)*hh(i,j,1)
      enddo
      enddo
      call mlwhalf(rosc,rosch,drosc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  energy ---
      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vy(i,j)**2
         ep=pr(i,j)*gm/(gm-1.)+0.5*ro(i,j)*vv
         fx(i,j)= ep*vx(i,j)-by(i,j)*ez(i,j)*pi4i
         fy(i,j)= ep*vy(i,j)+bx(i,j)*ez(i,j)*pi4i
         fx(i,j)= fx(i,j)*hh(i,j,2)*hh(i,j,3)
         fy(i,j)= fy(i,j)*hh(i,j,3)*hh(i,j,1)
      enddo
      enddo
      call mlwhalf(eesc,eesch,deesc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         t11=ro(i,j)*vx(i,j)**2+pr(i,j)
     &         -pi4i*bx(i,j)**2+pm(i,j)
         t22=ro(i,j)*vy(i,j)**2+pr(i,j)
     &         -pi4i*by(i,j)**2+pm(i,j)
         t33=pr(i,j)+pm(i,j)
         t12=ro(i,j)*vx(i,j)*vy(i,j)
     &         -pi4i*bx(i,j)*by(i,j)
         t21= t12
         t13=0.d0
         fx(i,j)= t11
         fy(i,j)= t21
         fx(i,j)= fx(i,j)*hh(i,j,2)*hh(i,j,3)
         fy(i,j)= fy(i,j)*hh(i,j,3)*hh(i,j,1)
         ss(i,j)=-gd(i,j,3)*t22-gi(i,j,2)*t33
     &           +gi(i,j,3)*t12+gd(i,j,2)*t13
         ss(i,j)= ss(i,j)*hh(i,j,1)*hh(i,j,2)*hh(i,j,3)
      enddo
      enddo
      call mlwhalf(rxsc,rxsch,drxsc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rxsch,drxsc,dt,ss,ix,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix
         t11=ro(i,j)*vx(i,j)**2+pr(i,j)
     &         -pi4i*bx(i,j)**2+pm(i,j)
         t22=ro(i,j)*vy(i,j)**2+pr(i,j)
     &         -pi4i*by(i,j)**2+pm(i,j)
         t33=pr(i,j)+pm(i,j)
         t23=0.d0
         t21=ro(i,j)*vy(i,j)*vx(i,j)
     &         -pi4i*by(i,j)*bx(i,j)
         t12=t21
         fx(i,j)= t12
         fy(i,j)= t22
         fx(i,j)= fx(i,j)*hh(i,j,2)*hh(i,j,3)
         fy(i,j)= fy(i,j)*hh(i,j,3)*hh(i,j,1)
         ss(i,j)=-gd(i,j,1)*t33-gi(i,j,3)*t11
     &           +gi(i,j,1)*t23+gd(i,j,3)*t21
         ss(i,j)= ss(i,j)*hh(i,j,1)*hh(i,j,2)*hh(i,j,3)
      enddo
      enddo
      call mlwhalf(rysc,rysch,drysc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rysch,drysc,dt,ss,ix,jx)

c---  x-magnetic ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= 0.
         fy(i,j)= ez(i,j)
         fy(i,j)= fy(i,j)*hh(i,j,3)
      enddo
      enddo
      call mlwhalf(bxsc,bxsch,dbxsc,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  y-magnetic ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= -ez(i,j)
         fy(i,j)= 0.
         fx(i,j)= fx(i,j)*hh(i,j,3)
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
         scm=hhm(i,j,1)*hhm(i,j,2)*hhm(i,j,3)
         roh(i,j) = rosch(i,j)/scm
         eeh      = eesch(i,j)/scm
         rxh      = rxsch(i,j)/scm
         ryh      = rysch(i,j)/scm
         vxh(i,j) = rxh/roh(i,j)
         vyh(i,j) = ryh/roh(i,j)
         bxh(i,j) = bxsch(i,j)/(hhm(i,j,2)*hhm(i,j,3))
         byh(i,j) = bysch(i,j)/(hhm(i,j,3)*hhm(i,j,1))
         ezh(i,j)=-vxh(i,j)*byh(i,j)+vyh(i,j)*bxh(i,j)
         bb=bxh(i,j)**2+byh(i,j)**2
         pmh(i,j) = pi8i*bb
         vv=vxh(i,j)**2+vyh(i,j)**2
         prh(i,j) = (gm-1)*(eeh-0.5*roh(i,j)*vv-pmh(i,j))
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
         fx(i,j)= fx(i,j)*hhm(i,j,2)*hhm(i,j,3)
         fy(i,j)= fy(i,j)*hhm(i,j,3)*hhm(i,j,1)
      enddo
      enddo
      call mlwfull(drosc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  energy     ---
      do j=1,jx-1
      do i=1,ix-1
         vv=vxh(i,j)**2+vyh(i,j)**2
         eph=prh(i,j)*gm/(gm-1.)+0.5*roh(i,j)*vv
         fx(i,j)= eph*vxh(i,j)-byh(i,j)*ezh(i,j)*pi4i
         fy(i,j)= eph*vyh(i,j)+bxh(i,j)*ezh(i,j)*pi4i
         fx(i,j)= fx(i,j)*hhm(i,j,2)*hhm(i,j,3)
         fy(i,j)= fy(i,j)*hhm(i,j,3)*hhm(i,j,1)
      enddo
      enddo
      call mlwfull(deesc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         t11=roh(i,j)*vxh(i,j)**2+prh(i,j)
     &          -pi4i*bxh(i,j)**2+pmh(i,j)
         t22=roh(i,j)*vyh(i,j)**2+prh(i,j)
     &          -pi4i*byh(i,j)**2+pmh(i,j)
         t33=prh(i,j)+pmh(i,j)
         t12=roh(i,j)*vxh(i,j)*vyh(i,j)
     &          -pi4i*bxh(i,j)*byh(i,j)
         t21=t12
         t13=0.d0
         fx(i,j)= t11
         fy(i,j)= t21
         fx(i,j)= fx(i,j)*hhm(i,j,2)*hhm(i,j,3)
         fy(i,j)= fy(i,j)*hhm(i,j,3)*hhm(i,j,1)
         ss(i,j)=-gdm(i,j,3)*t22-gim(i,j,2)*t33
     &           +gim(i,j,3)*t12+gdm(i,j,2)*t13
         ss(i,j)= ss(i,j)*hhm(i,j,1)*hhm(i,j,2)*hhm(i,j,3)
      enddo
      enddo
      call mlwfull(drxsc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drxsc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         t11=roh(i,j)*vxh(i,j)**2+prh(i,j)
     &          -pi4i*bxh(i,j)**2+pmh(i,j)
         t22=roh(i,j)*vyh(i,j)**2+prh(i,j)
     &          -pi4i*byh(i,j)**2+pmh(i,j)
         t33=prh(i,j)+pmh(i,j)
         t23=0.d0
         t21=roh(i,j)*vyh(i,j)*vxh(i,j)
     &          -pi4i*byh(i,j)*bxh(i,j)
         t12=t21
         fx(i,j)= t12
         fy(i,j)= t22
         fx(i,j)= fx(i,j)*hhm(i,j,2)*hhm(i,j,3)
         fy(i,j)= fy(i,j)*hhm(i,j,3)*hhm(i,j,1)
         ss(i,j)=-gdm(i,j,1)*t33-gim(i,j,3)*t11
     &           +gim(i,j,1)*t23+gdm(i,j,3)*t21
         ss(i,j)= ss(i,j)*hhm(i,j,1)*hhm(i,j,2)*hhm(i,j,3)
      enddo
      enddo
      call mlwfull(drysc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drysc,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  x-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= 0.
         fy(i,j)= ezh(i,j)
         fy(i,j)= fy(i,j)*hhm(i,j,3)
      enddo
      enddo
      call mlwfull(dbxsc,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  y-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= -ezh(i,j)
         fy(i,j)= 0.
         fx(i,j)= fx(i,j)*hhm(i,j,3)
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
      call mlwartv(eesc,deesc,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
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
         eesc(i,j) = eesc(i,j) +deesc(i,j)
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
         sc=hh(i,j,1)*hh(i,j,2)*hh(i,j,3)
         ro(i,j) = rosc(i,j)/sc
         ee      = eesc(i,j)/sc
         rx      = rxsc(i,j)/sc
         ry      = rysc(i,j)/sc
         vx(i,j) = rx/ro(i,j)
         vy(i,j) = ry/ro(i,j)
         bx(i,j) = bxsc(i,j)/(hh(i,j,2)*hh(i,j,3))
         by(i,j) = bysc(i,j)/(hh(i,j,3)*hh(i,j,1))
         vv=vx(i,j)**2+vy(i,j)**2
         bb=bx(i,j)**2+by(i,j)**2
         pr(i,j) = (gm-1)*(ee - 0.5*ro(i,j)*vv -pi8i*bb)
      enddo
      enddo

      return
      end
