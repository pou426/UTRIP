c======================================================================|
      subroutine mlw_m3t(ro,vx,vy,vz,by,bz,bx,bxm,dt,qav,cs2,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_m3t
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal 3-component MHD
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    vz(ix): [double] velocity 
c    by(ix): [double] magnetic field
c    bz(ix): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    cs2: [double] square of sound speed
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)

      dimension ro(ix),vx(ix),vy(ix),vz(ix),by(ix),bz(ix)
      dimension rx(ix),ry(ix),rz(ix)
      dimension ey(ix),ez(ix)

      dimension bx(ix),bxm(ix)

      dimension roh(ix),rxh(ix),ryh(ix),rzh(ix),byh(ix),bzh(ix)
      dimension vxh(ix),vyh(ix),vzh(ix)
      dimension eyh(ix),ezh(ix)

      dimension dro(ix),drx(ix),dry(ix),drz(ix),dby(ix),dbz(ix)

      dimension fx(ix),qx(ix)

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi8i=1./pi/8.
      pi4i=1./pi/4.
c----------------------------------------------------------------------|
c     ready
c----------------------------------------------------------------------|
      do i=1,ix
         dxim(i) = 1.0/dxm(i)
         dxi(i) = 1.0/dx(i)
      enddo

c----------------------------------------------------------------------|
c     initialize dro etc.                                   
c----------------------------------------------------------------------|

      do i=1,ix
         dro(i) = 0.0
         drx(i) = 0.0
         dry(i) = 0.0
         dby(i) = 0.0
         drz(i) = 0.0
         dbz(i) = 0.0
      enddo
c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         rx(i) = ro(i)*vx(i)
         ry(i) = ro(i)*vy(i)
         rz(i) = ro(i)*vz(i)
      enddo
      do i=1,ix
         ey(i) = -vz(i)*bx(i)+vx(i)*bz(i)
         ez(i) = -vx(i)*by(i)+vy(i)*bx(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)
      enddo
      call mlwhalf(ro ,roh ,dro ,dt,fx,dxi,dxim,ix)

c---  x-momentum ---
      do i=1,ix       
         bbm=by(i)**2+bz(i)**2
c        bbm=by(i)**2+bz(i)**2-bx(i)**2
         fx(i)= ro(i)*vx(i)**2+cs2*ro(i)+bbm*pi8i
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix)

c---  y-momentum ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)*vy(i)-bx(i)*by(i)*pi4i
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix)

c---  z-momentum ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)*vz(i)-bx(i)*bz(i)*pi4i
      enddo
      call mlwhalf(rz,rzh,drz,dt,fx,dxi,dxim,ix)

c---  y-magnetic ---
      do i=1,ix       
         fx(i)= -ez(i)
      enddo
      call mlwhalf(by ,byh ,dby ,dt,fx,dxi,dxim,ix)

c---  z-magnetic ---
      do i=1,ix       
         fx(i)= ey(i)
      enddo
      call mlwhalf(bz ,bzh ,dbz ,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         vxh(i)   = rxh(i)/roh(i)
         vyh(i)   = ryh(i)/roh(i)
         vzh(i)   = rzh(i)/roh(i)
      enddo
      do i=1,ix-1
         eyh(i) = -vzh(i)*bxm(i)+vxh(i)*bzh(i)
         ezh(i) = -vxh(i)*byh(i)+vyh(i)*bxm(i)
      enddo
c----------------------------------------------------------------------|
c     step intermediate results for full step                        
c----------------------------------------------------------------------|
c---  density  ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)
      enddo
      call mlwfull(dro ,dt,fx,dxi,ix)

c---  x-momentum ---
      do i=1,ix-1
         bbm=byh(i)**2+bzh(i)**2
c        bbm=byh(i)**2+bzh(i)**2-bxm(i)**2
         fx(i)= roh(i)*vxh(i)**2+cs2*roh(i)+bbm*pi8i
      enddo
      call mlwfull(drx,dt,fx,dxi,ix)

c---  y-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)*vyh(i)-bxm(i)*byh(i)*pi4i
      enddo
      call mlwfull(dry,dt,fx,dxi,ix)

c---  z-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)*vzh(i)-bxm(i)*bzh(i)*pi4i
      enddo
      call mlwfull(drz,dt,fx,dxi,ix)

c---  y-magnetic ---
      do i=1,ix-1
         fx(i)= -ezh(i)
      enddo
      call mlwfull(dby ,dt,fx,dxi,ix)

c---  z-magnetic ---
      do i=1,ix-1
         fx(i)=  eyh(i)
      enddo
      call mlwfull(dbz ,dt,fx,dxi,ix)
c----------------------------------------------------------------------|
c     diffusion coefficients for artificial viscosity             
c----------------------------------------------------------------------|
c     qav=3.0
      zero=0.0
      do i=1,ix-1
         qx(i)=qav*dxm(i)*max(zero,abs(vx(i+1)-vx(i))-1.0e-4)
      enddo
c----------------------------------------------------------------------|
c     apply artificial viscosity                                
c----------------------------------------------------------------------|
      call mlwartv(ro,dro,dt,qx,dxi,dxim,ix)
      call mlwartv(rx,drx,dt,qx,dxi,dxim,ix)
      call mlwartv(ry,dry,dt,qx,dxi,dxim,ix)
      call mlwartv(rz,drz,dt,qx,dxi,dxim,ix)
      call mlwartv(by,dby,dt,qx,dxi,dxim,ix)
      call mlwartv(bz,dbz,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = ro(i) +dro(i)
         rx(i) = rx(i) +drx(i)
         ry(i) = ry(i) +dry(i)
         rz(i) = rz(i) +drz(i)
         by(i) = by(i) +dby(i)
         bz(i) = bz(i) +dbz(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
         vy(i) = ry(i)/ro(i)
         vz(i) = rz(i)/ro(i)
      enddo

      return
      end
