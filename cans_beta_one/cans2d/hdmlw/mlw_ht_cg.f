c======================================================================|
      subroutine mlw_ht_cg(ro,vx,vz,dt,qav,cs2
     &         ,gx,gxm,gz,gzm,x,xm,dx,dxm,ix,dz,dzm,jx)
c======================================================================|
c
c NAME  mlw_h_cg
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * hydrodynamics
c        * Cylindrical coordinate, axis-symmetry
c        * gravity
c 
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vz(ix,jx): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c 
c    gx(ix,jx), gxm(ix,jx) : [double] gravity
c    gz(ix,jx), gzm(ix,jx) : [double] gravity
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
      dimension ro(ix,jx),vx(ix,jx),vz(ix,jx)
      dimension rx(ix,jx),rz(ix,jx)
      dimension roh(ix,jx),rxh(ix,jx),rzh(ix,jx)
      dimension vxh(ix,jx),vzh(ix,jx)
      dimension dro(ix,jx),drx(ix,jx),drz(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fz(ix,jx),qz(ix,jx)
      dimension ss(ix,jx)
      dimension gx(ix,jx), gxm(ix,jx)
      dimension gz(ix,jx), gzm(ix,jx)
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
         drz(i,j)= 0.0
      enddo
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         rx(i,j) = ro(i,j)*vx(i,j)
         rz(i,j) = ro(i,j)*vz(i,j)
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
         fx(i,j)= ro(i,j)*vx(i,j)**2+cs2*ro(i,j)
         fz(i,j)= ro(i,j)*vx(i,j)*vz(i,j)
         ss(i,j)= -ro(i,j)*vx(i,j)**2/x(i)+ro(i,j)*gx(i,j)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rxh ,drx ,dt,ss,ix,jx)

c---  z-momentum ---
      do j=1,jx
      do i=1,ix       
         fx(i,j)= ro(i,j)*vz(i,j)*vx(i,j)
         fz(i,j)= ro(i,j)*vz(i,j)**2+cs2*ro(i,j)
         ss(i,j)= -fx(i,j)/x(i)+ro(i,j)*gz(i,j)
      enddo
      enddo
      call mlwhalf(rz,rzh,drz,dt,fx,dxi,dxim,ix,fz,dzi,dzim,jx)
      call mlwsrch(rzh ,drz ,dt,ss,ix,jx)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         vxh(i,j)   = rxh(i,j)/roh(i,j)
         vzh(i,j)   = rzh(i,j)/roh(i,j)
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
         fx(i,j)= roh(i,j)*vxh(i,j)**2+cs2*roh(i,j)
         fz(i,j)= roh(i,j)*vxh(i,j)*vzh(i,j)
         ss(i,j)= -roh(i,j)*vxh(i,j)**2/xm(i)+roh(i,j)*gxm(i,j)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drx,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

c---  z-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vzh(i,j)*vxh(i,j)
         fz(i,j)= roh(i,j)*vzh(i,j)**2+cs2*roh(i,j)
         ss(i,j)= -fx(i,j)/xm(i)+roh(i,j)*gzm(i,j)
      enddo
      enddo
      call mlwfull(drz,dt,fx,dxi,ux0,ux1,ix,fz,dzi,uz0,uz1,jx)
      call mlwsrcf(drz,dt,ss,ux0,ux1,ix,uz0,uz1,jx)

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
      call mlwartv(rz,drz,dt,qx,dxi,dxim,ix,qz,dzi,dzim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         ro(i,j) = ro(i,j) +dro(i,j)
         rx(i,j) = rx(i,j) +drx(i,j)
         rz(i,j) = rz(i,j) +drz(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1  
         vx(i,j) = rx(i,j)/ro(i,j)
         vz(i,j) = rz(i,j)/ro(i,j)
      enddo
      enddo

      return
      end
