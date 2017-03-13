c======================================================================|
      subroutine mlw_h_tfix(ro,pr,vx,vy,dt,gm,gx,gxm,gy,gym
     &           ,x,dx,dxm,ix,y,dy,dym,jx,rkap,rkapm,visc,viscm
     &           ,hcs,hcsm,gasr,teb1,tebx)
c======================================================================|
c
c NAME  mlw_h_tfix
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * hydrodynamics
c        * gravity
c        * thermal conduction
c	 * viscosity
c        * y boundaries with fixed T and stress free 
c
c NOTE
c     margin must be 2 to use this module
c
c
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity 
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c 
c    gx(ix,jx), gxm(ix,jx) : [double] gravity
c    gy(ix,jx), gym(ix,jx) : [double] gravity
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c    rkap(ix,jx): [double] thermal conductivity
c    visc: [double] viscosity
c    gasr: [double] gas constant (pr=ro*te*gasr)
c    teb1: Boundary temperautre @ j=3
c    tebx: Boundary temperature @ j=jx-2 
c
c HISTORY
c    written 2002-6-21 H. Isobe based on mlh_h_g.f
c    
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)
      dimension x(ix),y(jx)
      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)
      dimension ux0(ix),ux1(ix)
      dimension dy(jx),dym(jx)
      dimension dyi(jx),dyim(jx)
      dimension uy0(jx),uy1(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension roh(ix,jx), rxh(ix,jx), eeh(ix,jx), ryh(ix,jx)
      dimension vxh(ix,jx), prh(ix,jx), vyh(ix,jx)
      dimension dro(ix,jx), dee(ix,jx), drx(ix,jx), dry(ix,jx)
      dimension fx(ix,jx),qx(ix,jx)
      dimension fy(ix,jx),qy(ix,jx)
      dimension ss(ix,jx)
      dimension gx(ix,jx), gxm(ix,jx)
      dimension gy(ix,jx), gym(ix,jx)
      dimension rkap(ix,jx),rkapm(ix,jx),visc(ix,jx),viscm(ix,jx)
      dimension te(ix,jx)
      dimension tdfx(ix,jx),tdfy(ix,jx)
      dimension visxx(ix,jx), visxy(ix,jx),visyy(ix,jx)
      dimension hcs(ix,jx),hcsm(ix,jx)

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
      do i=1,ix
         tdfx(i,j)=0.d0
         tdfy(i,j)=0.d0
         visxx(i,j)=0.d0
         visxy(i,j)=0.d0
         visyy(i,j)=0.d0
      enddo
      enddo
      

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         dro(i,j) = 0.0
         dee(i,j) = 0.0
         drx(i,j)= 0.0
         dry(i,j)= 0.0
      enddo
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         vv=vx(i,j)**2+vy(i,j)**2
         ee(i,j) = pr(i,j)/(gm-1)+0.5*ro(i,j)*vv
         rx(i,j) = ro(i,j)*vx(i,j)
         ry(i,j) = ro(i,j)*vy(i,j)
      enddo
      enddo

c----------------------------------------------------------------------|
c     viscous tensor
c----------------------------------------------------------------------|

      call visten(visxx,visxy,visyy,vx,vy,dx,ix,dy,jx,visc)

c----------------------------------------------------------------------|
c     thermal conduction flux
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         te(i,j)=pr(i,j)/ro(i,j)/gasr
      enddo
      enddo

      call cf_tfxh(tdfx,tdfy,te,dx,ix,dy,jx,rkap,teb1,tebx)


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

c---  energy ---

      do j=1,jx     
      do i=1,ix
         fx(i,j)= (ee(i,j)+pr(i,j))*vx(i,j) - tdfx(i,j)
     &           - visxx(i,j)*vx(i,j) - visxy(i,j)*vy(i,j)
         fy(i,j)= (ee(i,j)+pr(i,j))*vy(i,j) - tdfy(i,j)
     &           - visxy(i,j)*vx(i,j) - visyy(i,j)*vy(i,j)
         ss(i,j)= ro(i,j)*(vx(i,j)*gx(i,j)+vy(i,j)*gy(i,j))
     &            +hcs(i,j)         
      enddo
      enddo

      call mlwhalf(ee ,eeh ,dee ,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(eeh ,dee ,dt,ss,ix,jx)


c---  x-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vx(i,j)**2+pr(i,j)-visxx(i,j)
         fy(i,j)= ro(i,j)*vx(i,j)*vy(i,j)-visxy(i,j)
         ss(i,j)= ro(i,j)*gx(i,j)
      enddo
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)
      call mlwsrch(rxh ,drx ,dt,ss,ix,jx)

c---  y-momentum ---
      do j=1,jx
      do i=1,ix
         fx(i,j)= ro(i,j)*vy(i,j)*vx(i,j)-visxy(i,j)
         fy(i,j)= ro(i,j)*vy(i,j)**2+pr(i,j)-visyy(i,j)
         ss(i,j)= ro(i,j)*gy(i,j)
      enddo
      enddo

      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix,fy,dyi,dyim,jx)

      call mlwsrch(ryh ,dry ,dt,ss,ix,jx)

c----------------------------------------------------------------------|
c     calculate pressure from energy 
c----------------------------------------------------------------------|
      do j=1,jx-1
      do i=1,ix-1
         vxh(i,j)   = rxh(i,j)/roh(i,j)
         vyh(i,j)   = ryh(i,j)/roh(i,j)
         vv=vxh(i,j)**2+vyh(i,j)**2
         prh(i,j)   = (gm-1.d0)*(eeh(i,j)-0.5d0*roh(i,j)*vv)
      enddo
      enddo

c----------------------------------------------------------------------|
c     set the boundary conditions
c----------------------------------------------------------------------|

      do i=1,ix

         vyh(i,1)= 0.d0
         vyh(i,2)= 0.d0
         vxh(i,1)= vxh(i,4)
         vxh(i,2)= vxh(i,3)
         prh(i,1)= prh(i,4)
         prh(i,2)= prh(i,3)
         roh(i,1)= roh(i,4)
         roh(i,2)= roh(i,3)

	 vyh(i,jx-2)=0.d0
	 vyh(i,jx-1)=0.d0
	 vyh(i,jx  )=0.d0
         vxh(i,jx-2)= vxh(i,jx-3)
         vxh(i,jx-1)= vxh(i,jx-4)
         vxh(i,jx  )= vxh(i,jx-5)
         prh(i,jx-2)= prh(i,jx-3)
         prh(i,jx-1)= prh(i,jx-4)
         prh(i,jx  )= prh(i,jx-5)
         roh(i,jx-2)= roh(i,jx-3)
         roh(i,jx-1)= roh(i,jx-4)
         roh(i,jx  )= roh(i,jx-5)
      enddo

c----------------------------------------------------------------------|
c     viscous tensor for the 2nd step
c----------------------------------------------------------------------|

      call visten(visxx,visxy,visyy,vxh,vyh,dxim,ix,dyim,jx,viscm)

c----------------------------------------------------------------------|
c     thermal conduction flux for the 2nd step
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix-1
         te(i,j)=prh(i,j)/roh(i,j)/gasr
      enddo
      enddo

      call cf_tfxf(tdfx,tdfy,te,dxm,ix,dym,jx,rkapm,teb1,tebx)

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

c---  energy     ---

      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= (eeh(i,j)+prh(i,j))*vxh(i,j)-tdfx(i,j)
     &           -visxx(i,j)*vxh(i,j)-visxy(i,j)*vyh(i,j)
         fy(i,j)= (eeh(i,j)+prh(i,j))*vyh(i,j)-tdfy(i,j)
     &           -visxy(i,j)*vxh(i,j)-visyy(i,j)*vyh(i,j)
         ss(i,j)= roh(i,j)*(vxh(i,j)*gxm(i,j)+vyh(i,j)*gym(i,j))
     &            +hcsm(i,j)
      enddo
      enddo

      call mlwfull(dee ,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(dee,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  x-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vxh(i,j)**2+prh(i,j)-visxx(i,j)
         fy(i,j)= roh(i,j)*vxh(i,j)*vyh(i,j)-visxy(i,j)
         ss(i,j)= roh(i,j)*gxm(i,j)
      enddo
      enddo
      call mlwfull(drx,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(drx,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c---  y-momentum ---
      do j=1,jx-1
      do i=1,ix-1
         fx(i,j)= roh(i,j)*vyh(i,j)*vxh(i,j)-visxy(i,j)
         fy(i,j)= roh(i,j)*vyh(i,j)**2+prh(i,j)-visyy(i,j)
         ss(i,j)= roh(i,j)*gym(i,j)
      enddo
      enddo
      call mlwfull(dry,dt,fx,dxi,ux0,ux1,ix,fy,dyi,uy0,uy1,jx)
      call mlwsrcf(dry,dt,ss,ux0,ux1,ix,uy0,uy1,jx)

c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
      qav=1.0
      zero=0.0e0
      vcri=1.e-2

      do j=1,jx
      do i=1,ix-1
         qx(i,j)=qav*dxm(i)*max(zero,abs(vx(i+1,j)-vx(i,j))-vcri)
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
         qy(i,j)=qav*dym(j)*max(zero,abs(vy(i,j+1)-vy(i,j))-vcri)
      enddo
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(ro,dro,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(ee,dee,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix,qy,dyi,dyim,jx)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j) = ro(i,j) +dro(i,j)
         ee(i,j) = ee(i,j) +dee(i,j)
         rx(i,j)= rx(i,j)+drx(i,j)
         ry(i,j)= ry(i,j)+dry(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do j=3,jx-2  
      do i=3,ix-2  
         vx(i,j) = rx(i,j)/ro(i,j)
         vy(i,j) = ry(i,j)/ro(i,j)
         vv=vx(i,j)**2+vy(i,j)**2
         pr(i,j) = (gm-1)*(ee(i,j) - 0.5*ro(i,j)*vv)
      enddo
      enddo

      return
      end
