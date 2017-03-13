c======================================================================|
      subroutine mlw_m3t_c(ro,vx,vy,vz,bx,by,bz,ay,dt,qav,cs2
     &                              ,x,xm,dx,dxm,ix,dz,dzm,jx)
c======================================================================|
c
c NAME  mlw_m3t_c
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal 3-component MHD
c        * Cylindrical coordinate, axis-symmetry
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    vz(ix,jx): [double] velocity 
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    bz(ix,jx): [double] magnetic field
c    az(ix,jx): [double] magnetic vector potential
c    
c OUTPUTS
c    None
c 
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c    
c    dx(ix), dxm(ix): [double] grid spacing
c    dz(jx), dzm(jx): [double] grid spacing
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
      dimension x(ix),xm(ix)
      dimension dz(jx),dzm(jx)
      dimension dzi(jx),dzim(jx)
      dimension uz0(jx),uz1(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx),vz(ix,jx)
     &         ,bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension ex(ix,jx),ey(ix,jx),ez(ix,jx)
      dimension ay(ix,jx)
      dimension rx(ix,jx),ry(ix,jx),rz(ix,jx)
      dimension roh(ix,jx),rxh(ix,jx),ryh(ix,jx),rzh(ix,jx)
     &         ,bxh(ix,jx),byh(ix,jx),bzh(ix,jx)
      dimension exh(ix,jx),eyh(ix,jx),ezh(ix,jx)
      dimension ayh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx),vzh(ix,jx)
      dimension dro(ix,jx),drx(ix,jx),dry(ix,jx),drz(ix,jx)
     &         ,dbx(ix,jx),dby(ix,jx),dbz(ix,jx)
     &         ,day(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fz(ix,jx),qz(ix,jx)
      dimension ss(ix,jx)
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
         drx(i,j)= 0.0
         dry(i,j)= 0.0
         drz(i,j)= 0.0
         dbx(i,j)= 0.0
         dby(i,j)= 0.0
         dbz(i,j)= 0.0
         day(i,j)= 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
         rz(i,j) = ro(i,j)*vz(i,j)
         pr(i,j) = ro(i,j)*cs2
      enddo
      enddo
      do j=1,jx
      do i=1,ix
         ex(i,j) = -vy(i,j)*bz(i,j)+vz(i,j)*by(i,j)
         ey(i,j) = -vz(i,j)*bx(i,j)+vx(i,j)*bz(i,j)
         ez(i,j) = -vx(i,j)*by(i,j)+vy(i,j)*bx(i,j)
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

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)
     &          +pi8i*(by(i,j)**2+bz(i,j)**2-bx(i,j)**2)
         fz(i,j)= ro(i,j)*vx(i,j)*vz(i,j)-pi4i*bx(i,j)*bz(i,j)
         ss(i,j)= -(ro(i,j)*(vx(i,j)**2-vy(i,j)**2)
     &              +pi4i*(by(i,j)**2-bx(i,j)**2))/x(i)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rxh ,drx ,dt,ss,ix,jx)

c---  z-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= ro(i,j)*vz(i,j)*vx(i,j)-pi4i*bz(i,j)*bx(i,j)
         fz(i,j)= ro(i,j)*vz(i,j)**2+pr(i,j)
     &          +pi8i*(bx(i,j)**2+by(i,j)**2-bz(i,j)**2)
         ss(i,j)= -fx(i,j)/x(i)
      enddo
      enddo
      call mlwhalf(rz,rzh,drz,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rzh ,drz ,dt,ss,ix,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)-pi4i*by(i,j)*bx(i,j)
         fz(i,j)= ro(i,j)*vz(i,j)*vy(i,j)-pi4i*bz(i,j)*by(i,j)
         ss(i,j)= -fx(i,j)*2/x(i)
      enddo
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(ryh ,dry ,dt,ss,ix,jx)

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

c---  y-magnetic ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= -ez(i,j)
         fz(i,j)=  ex(i,j)
      enddo
      enddo
      call mlwhalf(by,byh,dby,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)

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
         vyh(i,j)   = ryh(i,j)/roh(i,j)
         vzh(i,j)   = rzh(i,j)/roh(i,j)
         prh(i,j)   = roh(i,j)*cs2
      enddo
      enddo
      do j=1,jx
      do i=1,ix
        exh(i,j)=-vyh(i,j)*bzh(i,j)+vzh(i,j)*byh(i,j)
        eyh(i,j)=-vzh(i,j)*bxh(i,j)+vxh(i,j)*bzh(i,j)
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
         fz(i,j)= roh(i,j)*vzh(i,j)
         ss(i,j)= -fx(i,j)/xm(i)
      enddo
      enddo
      call mlwfull(dro ,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(dro,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+prh(i,j)
     &          +pi8i*(byh(i,j)**2+bzh(i,j)**2-bxh(i,j)**2)
         fz(i,j)= roh(i,j)*vxh(i,j)*vzh(i,j)-pi4i*bxh(i,j)*bzh(i,j)
         ss(i,j)= -(roh(i,j)*(vxh(i,j)**2-vyh(i,j)**2)
     &              +pi4i*(byh(i,j)**2-bxh(i,j)**2))/xm(i)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drx,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  z-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vzh(i,j)*vxh(i,j)-pi4i*bzh(i,j)*bxh(i,j)
         fz(i,j)= roh(i,j)*vzh(i,j)**2+prh(i,j)
     &          +pi8i*(bxh(i,j)**2+byh(i,j)**2-bzh(i,j)**2)
         ss(i,j)= -fx(i,j)/xm(i)
      enddo
      enddo
      call mlwfull(drz,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drz,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= roh(i,j)*vyh(i,j)*vxh(i,j)-pi4i*byh(i,j)*bxh(i,j)
         fz(i,j)= roh(i,j)*vzh(i,j)*vyh(i,j)-pi4i*bzh(i,j)*byh(i,j)
         ss(i,j)= -fx(i,j)*2/xm(i)
      enddo
      enddo
      call mlwfull(dry,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(dry,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

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

c---  y-magnetic ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= -ezh(i,j)
         fz(i,j)=  exh(i,j)
      enddo
      enddo
      call mlwfull(dby,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)

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
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(rz,drz,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(bx,dbx,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(by,dby,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
      call mlwartv(bz,dbz,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         ro(i,j) = ro(i,j) +dro(i,j)
         rx(i,j) = rx(i,j) +drx(i,j)
         ry(i,j) = ry(i,j) +dry(i,j)
         rz(i,j) = rz(i,j) +drz(i,j)
         bx(i,j) = bx(i,j) +dbx(i,j)
         by(i,j) = by(i,j) +dby(i,j)
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
         vy(i,j) = ry(i,j)/ro(i,j)
         vz(i,j) = rz(i,j)/ro(i,j)
      enddo
      enddo

      return
      end
