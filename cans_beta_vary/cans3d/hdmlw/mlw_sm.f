c======================================================================|
      subroutine mlw_sm(ro,pr,vx,vy,vz,bx,by,bz,dt,qav,gm
     &                 ,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dxm(ix),dxi(ix),dxim(ix),ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx),dyi(jx),dyim(jx),uy0(jx),uy1(jx)
      dimension dz(kx),dzm(kx),dzi(kx),dzim(kx),uz0(kx),uz1(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx),de(ix,jx,kx),ee(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension roh(ix,jx,kx),prh(ix,jx,kx),deh(ix,jx,kx),eeh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension rxh(ix,jx,kx),ryh(ix,jx,kx),rzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)
      dimension exh(ix,jx,kx),eyh(ix,jx,kx),ezh(ix,jx,kx)
      dimension dde(ix,jx,kx),dee(ix,jx,kx)
      dimension drx(ix,jx,kx),dry(ix,jx,kx),drz(ix,jx,kx)
      dimension dbx(ix,jx,kx),dby(ix,jx,kx),dbz(ix,jx,kx)
      dimension fx(ix,jx,kx),qx(ix,jx,kx)
      dimension fy(ix,jx,kx),qy(ix,jx,kx)
      dimension fz(ix,jx,kx),qz(ix,jx,kx)

      dimension gl(ix,jx,kx),el(ix,jx,kx),pm(ix,jx,kx)
      dimension glm(ix,jx,kx),elh(ix,jx,kx),pmh(ix,jx,kx)

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
         dde(i,j,k)= 0.d0
         dee(i,j,k)= 0.d0
         drx(i,j,k)= 0.d0
         dry(i,j,k)= 0.d0
         drz(i,j,k)= 0.d0
         dbx(i,j,k)= 0.d0
         dby(i,j,k)= 0.d0
         dbz(i,j,k)= 0.d0
         glm(i,j,k)= 0.d0
         roh(i,j,k)= 0.d0
         prh(i,j,k)= 0.d0
         vxh(i,j,k)= 0.d0
         vyh(i,j,k)= 0.d0
         vzh(i,j,k)= 0.d0
         bxh(i,j,k)= 0.d0
         byh(i,j,k)= 0.d0
         bzh(i,j,k)= 0.d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         vf2=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         gl(i,j,k)=1.d0/sqrt(1.d0-vf2)
         ex(i,j,k) = -vy(i,j,k)*bz(i,j,k)+vz(i,j,k)*by(i,j,k)
         ey(i,j,k) = -vz(i,j,k)*bx(i,j,k)+vx(i,j,k)*bz(i,j,k)
         ez(i,j,k) = -vx(i,j,k)*by(i,j,k)+vy(i,j,k)*bx(i,j,k)
         ef2=ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2
         bf2=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         pm(i,j,k) = pi8i*(bf2+ef2)
         de(i,j,k) = gl(i,j,k)*ro(i,j,k)
         en = ro(i,j,k)+pr(i,j,k)/(gm-1)
         el(i,j,k) = (en+pr(i,j,k))*gl(i,j,k)**2
         ee(i,j,k) = el(i,j,k)-pr(i,j,k)-de(i,j,k)
         rx(i,j,k) = el(i,j,k)*vx(i,j,k)
         ry(i,j,k) = el(i,j,k)*vy(i,j,k)
         rz(i,j,k) = el(i,j,k)*vz(i,j,k)
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
         fx(i,j,k)= de(i,j,k)*vx(i,j,k)
         fy(i,j,k)= de(i,j,k)*vy(i,j,k)
         fz(i,j,k)= de(i,j,k)*vz(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(de,deh,dde,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  energy ---
      do k=1,kx
      do j=1,jx       
      do i=1,ix       
         fx(i,j,k)= rx(i,j,k)-de(i,j,k)*vx(i,j,k)
         fy(i,j,k)= ry(i,j,k)-de(i,j,k)*vy(i,j,k)
         fz(i,j,k)= rz(i,j,k)-de(i,j,k)*vz(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(ee,eeh,dee,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  x-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= el(i,j,k)*vx(i,j,k)**2+pr(i,j,k)
     &                  -pi4i*bx(i,j,k)**2+pm(i,j,k)
         fy(i,j,k)= el(i,j,k)*vx(i,j,k)*vy(i,j,k)
     &                  -pi4i*bx(i,j,k)*by(i,j,k)
         fz(i,j,k)= el(i,j,k)*vx(i,j,k)*vz(i,j,k)
     &                  -pi4i*bx(i,j,k)*bz(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  y-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= el(i,j,k)*vy(i,j,k)*vx(i,j,k)
     &                  -pi4i*by(i,j,k)*bx(i,j,k)
         fy(i,j,k)= el(i,j,k)*vy(i,j,k)**2+pr(i,j,k)
     &                  -pi4i*by(i,j,k)**2+pm(i,j,k)
         fz(i,j,k)= el(i,j,k)*vy(i,j,k)*vz(i,j,k)
     &                  -pi4i*by(i,j,k)*bz(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(ry,ryh,dry,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  z-momentum ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)= el(i,j,k)*vz(i,j,k)*vx(i,j,k)
     &                  -pi4i*bz(i,j,k)*bx(i,j,k)
         fy(i,j,k)= el(i,j,k)*vz(i,j,k)*vy(i,j,k)
     &                  -pi4i*bz(i,j,k)*by(i,j,k)
         fz(i,j,k)= el(i,j,k)*vz(i,j,k)**2+pr(i,j,k)
     &                  -pi4i*bz(i,j,k)**2+pm(i,j,k)
      enddo
      enddo
      enddo
      call mlwhalf(rz,rzh,drz,dt
     &             ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c---  x-magnetic ---
      do k=1,kx
      do j=1,jx
      do i=1,ix
         fx(i,j,k)=  0.
         fy(i,j,k)=  ez(i,j,k)
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
         fy(i,j,k)= -ex(i,j,k)
         fz(i,j,k)=  0.
      enddo
      enddo
      enddo
      call mlwhalf(bz,bzh,dbz,dt
     &      ,fx,dxi,dxim,ix,fy,dyi,dyim,jx,fz,dzi,dzim,kx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      mix=100
      tolf=1.d-20
      tolx=1.d-20
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         roh(i,j,k) = ro(i,j,k)
         prh(i,j,k) = pr(i,j,k)
         vxh(i,j,k) = vx(i,j,k)
         vyh(i,j,k) = vy(i,j,k)
         vzh(i,j,k) = vz(i,j,k)
      enddo
      enddo
      enddo

      call rtnewt_sm(roh,prh,vxh,vyh,vzh,glm
     &     ,deh,eeh,rxh,ryh,rzh,bxh,byh,bzh,gm,mix,tolf,tolx
     &     ,ix,1,ix-1,jx,1,jx-1,kx,1,kx-1)

      do k=1,kx
      do j=1,jx
      do i=1,ix
         exh(i,j,k) = -vyh(i,j,k)*bzh(i,j,k)+vzh(i,j,k)*byh(i,j,k)
         eyh(i,j,k) = -vzh(i,j,k)*bxh(i,j,k)+vxh(i,j,k)*bzh(i,j,k)
         ezh(i,j,k) = -vxh(i,j,k)*byh(i,j,k)+vyh(i,j,k)*bxh(i,j,k)
         ef2=exh(i,j,k)**2+eyh(i,j,k)**2+ezh(i,j,k)**2
         bf2=bxh(i,j,k)**2+byh(i,j,k)**2+bzh(i,j,k)**2
         pmh(i,j,k) = pi8i*(bf2+ef2)
         enh = roh(i,j,k)+prh(i,j,k)/(gm-1)
         elh(i,j,k) = (enh+prh(i,j,k))*glm(i,j,k)**2
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
         fx(i,j,k)= deh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= deh(i,j,k)*vyh(i,j,k)
         fz(i,j,k)= deh(i,j,k)*vzh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dde,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  energy     ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= rxh(i,j,k)-deh(i,j,k)*vxh(i,j,k)
         fy(i,j,k)= ryh(i,j,k)-deh(i,j,k)*vyh(i,j,k)
         fz(i,j,k)= rzh(i,j,k)-deh(i,j,k)*vzh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dee,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  x-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= elh(i,j,k)*vxh(i,j,k)**2+prh(i,j,k)
     &                  -pi4i*bxh(i,j,k)**2+pmh(i,j,k)
         fy(i,j,k)= elh(i,j,k)*vxh(i,j,k)*vyh(i,j,k)
     &                  -pi4i*bxh(i,j,k)*byh(i,j,k)
         fz(i,j,k)= elh(i,j,k)*vxh(i,j,k)*vzh(i,j,k)
     &                  -pi4i*bxh(i,j,k)*bzh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(drx,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  y-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= elh(i,j,k)*vyh(i,j,k)*vxh(i,j,k)
     &                  -pi4i*byh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= elh(i,j,k)*vyh(i,j,k)**2+prh(i,j,k)
     &                  -pi4i*byh(i,j,k)**2+pmh(i,j,k)
         fz(i,j,k)= elh(i,j,k)*vyh(i,j,k)*vzh(i,j,k)
     &                  -pi4i*byh(i,j,k)*bzh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(dry,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  z-momentum ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)= elh(i,j,k)*vzh(i,j,k)*vxh(i,j,k)
     &                  -pi4i*bzh(i,j,k)*bxh(i,j,k)
         fy(i,j,k)= elh(i,j,k)*vzh(i,j,k)*vyh(i,j,k)
     &                  -pi4i*bzh(i,j,k)*byh(i,j,k)
         fz(i,j,k)= elh(i,j,k)*vzh(i,j,k)**2+prh(i,j,k)
     &                  -pi4i*bzh(i,j,k)**2+pmh(i,j,k)
      enddo
      enddo
      enddo
      call mlwfull(drz,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c---  x-magnetic ---
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j,k)=  0.
         fy(i,j,k)=  ezh(i,j,k)
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
         fy(i,j,k)= -exh(i,j,k)
         fz(i,j,k)=  0.
      enddo
      enddo
      enddo
      call mlwfull(dbz,dt
     &       ,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx,fz,dzi,uz0,uz1,kx)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=1.0
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
      call mlwartv(de,dde,dt
     &              ,qx,dxi,dxim,ix,qy,dyi,dyim,jx,qz,dzi,dzim,kx)
      call mlwartv(ee,dee,dt
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
         de(i,j,k)= de(i,j,k)+dde(i,j,k)
         ee(i,j,k)= ee(i,j,k)+dee(i,j,k)
         rx(i,j,k)= rx(i,j,k)+drx(i,j,k)
         ry(i,j,k)= ry(i,j,k)+dry(i,j,k)
         rz(i,j,k)= rz(i,j,k)+drz(i,j,k)
         bx(i,j,k)= bx(i,j,k)+dbx(i,j,k)
         by(i,j,k)= by(i,j,k)+dby(i,j,k)
         bz(i,j,k)= bz(i,j,k)+dbz(i,j,k)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      call rtnewt_sm(ro,pr,vx,vy,vz,gl
     &     ,de,ee,rx,ry,rz,bx,by,bz,gm,mix,tolf,tolx
     &     ,ix,2,ix-1,jx,2,jx-1,kx,2,kx-1)


      return
      end
