c======================================================================|
      subroutine mlw_ht(ro,vx,vy,dt,qav,cs2,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  mlw_ht
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal hydrodynamics
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity 
c    
c OUTPUTS
c    None
c 
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c 
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
      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension rx(ix,jx),ry(ix,jx)
      dimension roh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension vxh(ix,jx),vyh(ix,jx)
      dimension dro(ix,jx),drx(ix,jx),dry(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)
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
         dro(i,j) = 0.0
         drx(i,j)= 0.0
         dry(i,j)= 0.0
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
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
      enddo
      enddo
      call mlwhalf(ro ,roh ,dro,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+cs2*ro(i,j)
         fy(i,j)= ro(i,j)*vx(i,j)*vy(i,j)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)**2+cs2*ro(i,j)
      enddo
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         vxh(i,j)   = rxh(i,j)/roh(i,j)
         vyh(i,j)   = ryh(i,j)/roh(i,j)
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
      enddo
      enddo
      call mlwfull(dro ,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+cs2*roh(i,j)
         fy(i,j)= roh(i,j)*vxh(i,j)*vyh(i,j)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vyh(i,j)*vxh(i,j)
         fy(i,j)= roh(i,j)*vyh(i,j)**2+cs2*roh(i,j)
      enddo
      enddo
      call mlwfull(dry,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
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
      call mlwartv(ro,dro,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         ro(i,j) = ro(i,j) +dro(i,j)
         rx(i,j) = rx(i,j) +drx(i,j)
         ry(i,j) = ry(i,j) +dry(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=2,jx-1  
      do i=2,ix-1  
         vx(i,j) = rx(i,j)/ro(i,j)
         vy(i,j) = ry(i,j)/ro(i,j)
      enddo
      enddo

      return
      end
