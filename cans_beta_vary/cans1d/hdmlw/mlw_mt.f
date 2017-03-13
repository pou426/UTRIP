c======================================================================|
      subroutine mlw_mt(ro,vx,vy,by,bx,bxm,dt,qav,cs2,dx,dxm,ix)
c======================================================================|
c
c NAME  mlw_mt
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * isothermal MHD
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vx(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    by(ix): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    dx(ix), dxm(ix): [double] grid spacing
c    cs2: [double] square of sound speed
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dxi(ix),dxim(ix)

      dimension ro(ix),vx(ix),vy(ix),by(ix)
      dimension rx(ix),ry(ix)
      dimension ez(ix)

      dimension bx(ix),bxm(ix)

      dimension roh(ix),rxh(ix),ryh(ix),byh(ix)
      dimension vxh(ix),vyh(ix)
      dimension ezh(ix)

      dimension dro(ix),drx(ix),dry(ix),dby(ix)

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
      enddo

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      do i=1,ix
         rx(i) = ro(i)*vx(i)
         ry(i) = ro(i)*vy(i)
      enddo
      do i=1,ix
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
         bbm=by(i)**2
c        bbm=by(i)**2-bx(i)**2
         fx(i)= ro(i)*vx(i)**2+cs2*ro(i)+bbm*pi8i
      enddo
      call mlwhalf(rx,rxh,drx,dt,fx,dxi,dxim,ix)

c---  y-momentum ---
      do i=1,ix       
         fx(i)= ro(i)*vx(i)*vy(i)-bx(i)*by(i)*pi4i
      enddo
      call mlwhalf(ry,ryh,dry,dt,fx,dxi,dxim,ix)

c---  y-magnetic ---
      do i=1,ix       
         fx(i)= -ez(i)
      enddo
      call mlwhalf(by ,byh ,dby ,dt,fx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=1,ix-1
         vxh(i)   = rxh(i)/roh(i)
         vyh(i)   = ryh(i)/roh(i)
      enddo
      do i=1,ix
        ezh(i)=-vxh(i)*byh(i)+vyh(i)*bxm(i)
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
         bbm=byh(i)**2
c        bbm=byh(i)**2-bxm(i)**2
         fx(i)= roh(i)*vxh(i)**2+cs2*roh(i)+bbm*pi8i
      enddo
      call mlwfull(drx,dt,fx,dxi,ix)

c---  y-momentum ---
      do i=1,ix-1
         fx(i)= roh(i)*vxh(i)*vyh(i)-bxm(i)*byh(i)*pi4i
      enddo
      call mlwfull(dry,dt,fx,dxi,ix)

c---  y-magnetic ---
      do i=1,ix-1
         fx(i)= -ezh(i)
      enddo
      call mlwfull(dby ,dt,fx,dxi,ix)
c-------------------------------------------------------------------|
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
      call mlwartv(by,dby,dt,qx,dxi,dxim,ix)
c----------------------------------------------------------------------|
c     update internal points                           
c----------------------------------------------------------------------|
      do i=2,ix-1
         ro(i) = ro(i) +dro(i)
         rx(i) = rx(i) +drx(i)
         ry(i) = ry(i) +dry(i)
         by(i) = by(i) +dby(i)
      enddo
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
      do i=2,ix-1  
         vx(i) = rx(i)/ro(i)
         vy(i) = ry(i)/ro(i)
      enddo

      return
      end
