c======================================================================|
      subroutine mlw_mt_cg(ro,vx,vy,vz,bx,by,bz,dt,qav,cs2
     &                 ,gx,gxm,gy,gym,gz,gzm
     &                 ,x,xm,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension roh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)
      dimension exh(ix,jx,kx),eyh(ix,jx,kx),ezh(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension rxh(ix,jx,kx),ryh(ix,jx,kx),rzh(ix,jx,kx)
      dimension dro(ix,jx,kx)
      dimension drx(ix,jx,kx),dry(ix,jx,kx),drz(ix,jx,kx)
      dimension dbx(ix,jx,kx),dby(ix,jx,kx),dbz(ix,jx,kx)
      dimension fx(ix,jx,kx),qx(ix,jx,kx)
      dimension fy(ix,jx,kx),qy(ix,jx,kx)
      dimension fz(ix,jx,kx),qz(ix,jx,kx)
      dimension ss(ix,jx,kx)
      dimension x(ix),xm(ix)
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

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         dro(i,j,k)= 0.0
         drx(i,j,k)= 0.0
         dry(i,j,k)= 0.0
         drz(i,j,k)= 0.0
         dbx(i,j,k)= 0.0
         dby(i,j,k)= 0.0
         dbz(i,j,k)= 0.0
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
         ex(i,j,k) = -vy(i,j,k)*bz(i,j,k)+vz(i,j,k)*by(i,j,k)
         ey(i,j,k) = -vz(i,j,k)*bx(i,j,k)+vx(i,j,k)*bz(i,j,k)
         ez(i,j,k) = -vx(i,j,k)*by(i,j,k)+vy(i,j,k)*bx(i,j,k)
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
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)/x(i)
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)
         ss(i,j,k)= -fx(i,j,k)/x(i)
      enddo
      enddo
      enddo
      call mlwhalf(ro,roh,dro,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(roh,dro,dt,ss,ix,jx,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vx(i,j,k)**2+cs2*ro(i,j,k)
     &          +pi8i*(by(i,j,k)**2+bz(i,j,k)**2-bx(i,j,k)**2)
         fy(i,j,k)= ro(i,j,k)*vx(i,j,k)*vy(i,j,k)
     &          -pi4i*bx(i,j,k)*by(i,j,k)
         fy(i,j,k)= fy(i,j,k)/x(i)
         fz(i,j,k)= ro(i,j,k)*vx(i,j,k)*vz(i,j,k)
     &          -pi4i*bx(i,j,k)*bz(i,j,k)
         ss(i,j,k)= -ro(i,j,k)/x(i)*(vx(i,j,k)**2-vy(i,j,k)**2)
     &                   +pi4i/x(i)*(bx(i,j,k)**2-by(i,j,k)**2)
     &             +ro(i,j,k)*gx(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rxh,drx,dt,ss,ix,jx,kx)

c---  y-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vy(i,j,k)*vx(i,j,k)
     &          -pi4i*by(i,j,k)*bx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vy(i,j,k)**2+cs2*ro(i,j,k)
     &          +pi8i*(bx(i,j,k)**2+bz(i,j,k)**2-by(i,j,k)**2)
         fy(i,j,k)= fy(i,j,k)/x(i)
         fz(i,j,k)= ro(i,j,k)*vy(i,j,k)*vz(i,j,k)
     &          -pi4i*by(i,j,k)*bz(i,j,k)
         ss(i,j,k)= -2.d0*ro(i,j,k)/x(i)*vx(i,j,k)*vy(i,j,k)
     &                   +2.d0*pi4i/x(i)*bx(i,j,k)*by(i,j,k)
     &             +ro(i,j,k)*gy(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(ry,ryh,dry,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(ryh,dry,dt,ss,ix,jx,kx)

c---  z-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= ro(i,j,k)*vz(i,j,k)*vx(i,j,k)
     &          -pi4i*bz(i,j,k)*bx(i,j,k)
         fy(i,j,k)= ro(i,j,k)*vz(i,j,k)*vy(i,j,k)
     &          -pi4i*bz(i,j,k)*by(i,j,k)
         fy(i,j,k)= fy(i,j,k)/x(i)
         fz(i,j,k)= ro(i,j,k)*vz(i,j,k)**2+cs2*ro(i,j,k)
     &          +pi8i*(bx(i,j,k)**2+by(i,j,k)**2-bz(i,j,k)**2)
         ss(i,j,k)= -fx(i,j,k)/x(i)
     &             +ro(i,j,k)*gz(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(rz,rzh,drz,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(rzh,drz,dt,ss,ix,jx,kx)

c---  x-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  0.
         fy(i,j,k)=  ez(i,j,k)/x(i)
         fz(i,j,k)= -ey(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(bx,bxh,dbx,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  y-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= -ez(i,j,k)
         fy(i,j,k)=  0.
         fz(i,j,k)=  ex(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(by,byh,dby,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  z-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  ey(i,j,k)
         fy(i,j,k)= -ex(i,j,k)/x(i)
         fz(i,j,k)=  0.
         ss(i,j,k)=  -fx(i,j,k)/x(i)
      enddo
      enddo
      enddo
      call mlwhalf(bz,bzh,dbz,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)
      call mlwsrch(bzh,dbz,dt,ss,ix,jx,kx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         vxh(i,j,k)   = rxh(i,j,k)/roh(i,j,k)
         vyh(i,j,k)   = ryh(i,j,k)/roh(i,j,k)
         vzh(i,j,k)   = rzh(i,j,k)/roh(i,j,k)
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
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= roh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vyh(i,j,k)/xm(i)
         fz(i,j,k)= roh(i,j,k)*vzh(i,j,k)
         ss(i,j,k)= -fx(i,j,k)/xm(i)
      enddo
      enddo
      enddo
      call mlwfull(dro,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(dro,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  x-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= roh(i,j,k)*vxh(i,j,k)**2+cs2*roh(i,j,k)
     &          +pi8i*(byh(i,j,k)**2+bzh(i,j,k)**2-bxh(i,j,k)**2)
         fy(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
     &          -pi4i*bxh(i,j,k)*byh(i,j,k)
         fy(i,j,k)= fy(i,j,k)/xm(i)
         fz(i,j,k)= roh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
     &          -pi4i*bxh(i,j,k)*bzh(i,j,k)
         ss(i,j,k)= -roh(i,j,k)/xm(i)*(vxh(i,j,k)**2-vyh(i,j,k)**2)
     &                    +pi4i/xm(i)*(bxh(i,j,k)**2-byh(i,j,k)**2)
     &             +roh(i,j,k)*gxm(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(drx,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drx,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  y-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= roh(i,j,k)*vyh(i,j,k)*vxh(i,j,k)
     &          -pi4i*byh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vyh(i,j,k)**2+cs2*roh(i,j,k)
     &          +pi8i*(bxh(i,j,k)**2+bzh(i,j,k)**2-byh(i,j,k)**2)
         fy(i,j,k)= fy(i,j,k)/xm(i)
         fz(i,j,k)= roh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
     &          -pi4i*byh(i,j,k)*bzh(i,j,k)
         ss(i,j,k)= -2.d0*roh(i,j,k)/xm(i)*vxh(i,j,k)*vyh(i,j,k)
     &                    +2.d0*pi4i/xm(i)*bxh(i,j,k)*byh(i,j,k)
     &             +roh(i,j,k)*gym(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dry,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(dry,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  z-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= roh(i,j,k)*vzh(i,j,k)*vxh(i,j,k)
     &          -pi4i*bzh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= roh(i,j,k)*vzh(i,j,k)*vyh(i,j,k)
     &          -pi4i*bzh(i,j,k)*byh(i,j,k)
         fy(i,j,k)= fy(i,j,k)/xm(i)
         fz(i,j,k)= roh(i,j,k)*vzh(i,j,k)**2+cs2*roh(i,j,k)
     &          +pi8i*(bxh(i,j,k)**2+byh(i,j,k)**2-bzh(i,j,k)**2)
         ss(i,j,k)= -fx(i,j,k)/xm(i)
     &             +roh(i,j,k)*gzm(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(drz,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(drz,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c---  x-magnetic ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)=  0.
         fy(i,j,k)=  ezh(i,j,k)/xm(i)
         fz(i,j,k)= -eyh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dbx,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  y-magnetic ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= -ezh(i,j,k)
         fy(i,j,k)=  0.
         fz(i,j,k)=  exh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dby,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  z-magnetic ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)=  eyh(i,j,k)
         fy(i,j,k)= -exh(i,j,k)/xm(i)
         fz(i,j,k)=  0.
         ss(i,j,k)= -fx(i,j,k)/xm(i)
      enddo
      enddo
      enddo
      call mlwfull(dbz,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)
      call mlwsrcf(dbz,dt,ss,ux0,ux1,ix,uy0,uy1,jx,uz0,uz1,kx)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0e0
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
      call mlwartv(ro,dro,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rx,drx,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(ry,dry,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(rz,drz,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(bx,dbx,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(by,dby,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(bz,dbz,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         ro(i,j,k)= ro(i,j,k) +dro(i,j,k)
         rx(i,j,k)= rx(i,j,k) +drx(i,j,k)
         ry(i,j,k)= ry(i,j,k) +dry(i,j,k)
         rz(i,j,k)= rz(i,j,k) +drz(i,j,k)
         bx(i,j,k)= bx(i,j,k) +dbx(i,j,k)
         by(i,j,k)= by(i,j,k) +dby(i,j,k)
         bz(i,j,k)= bz(i,j,k) +dbz(i,j,k)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do k=2,kx-1
      do j=2,jx-1  
      do i=2,ix-1  
         vx(i,j,k) = rx(i,j,k)/ro(i,j,k)
         vy(i,j,k) = ry(i,j,k)/ro(i,j,k)
         vz(i,j,k) = rz(i,j,k)/ro(i,j,k)
      enddo
      enddo
      enddo

      return
      end
