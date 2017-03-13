c======================================================================|
      subroutine cip_m(ro,pr,vx,vy,by,te,vxm,rodx,tedx,vxdxm,vydx
     &    ,bx,bxm,dt,cvis,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_m
c
c PURPOSE
c    solve eqs. by CIP-MOCCT  method with effects of
c        * MHD
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    te(ix): [double] temperature
c    vxm(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    by(ix): [double] magnetic field
c    rodx(ix): [double] density gradient
c    tedx(ix): [double] temperature gradient
c    vxdxm(ix): [double] velocity gradient
c    vydx(ix): [double] velocity gradient
c
c OUTPUTS
c    vx(ix): [double] velocity 
c    pr(ix): [double] pressure
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    dx(ix), dxm(ix) : [double] grid spacing
c    gm: [double] polytropic index gamma
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)

      dimension  ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
      dimension  te(ix),vxm(ix)
      dimension  bx(ix),bxm(ix)
      dimension  roh(ix),teh(ix),vxmh(ix),vyh(ix),byh(ix)
      dimension  rodx(ix),tedx(ix),vxdxm(ix),vydx(ix)
      dimension  rodxh(ix),tedxh(ix),vxdxmh(ix),vydxh(ix)
      dimension  vyc(ix),byc(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c---  calculate specific thermal energy from pressure

      do i=1,ix
        roh(i) = ro(i)
        teh(i) = te(i)
        vxmh(i) = vxm(i)
        vyh(i) = vy(i)
        byh(i) = by(i)
        rodxh(i)=rodx(i)
        tedxh(i)=tedx(i)
        vxdxmh(i)=vxdxm(i)
        vydxh(i)=vydx(i)
      enddo

c     calculate pressure
      do i=1,ix
         pr(i) = ro(i)*te(i)/gm
      enddo
c--- velocity on the grid points ---
      do i=2,ix
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|

c     cvis=0.75
      do i=2,ix
        du=vxm(i)-vxm(i-1)
        if (du.lt.0.0) then
          cs=sqrt(te(i))
          pr(i)=pr(i)+cvis*(-ro(i)*cs*du+(gm+1)/2*ro(i)*du**2)
        endif
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### Calculate magnetc stress term Using MOC method

      call moclag(byc,ro,bxm,vy,by,dt,dxm,ix)
      do i=2,ix-1
        vyh(i)=vyh(i)+1./(4.*pi*ro(i))*bx(i)*(byc(i)-byc(i-1))/dx(i)*dt
      enddo

c--- ### important ### vx must be advanced first.
      do i=1,ix-1
        vxmh(i)=vxmh(i)+dt*(
     &    -((pr(i+1)-pr(i)))
     &        /dxm(i)/((ro(i)+ro(i+1))/2)
     &    -(pi4i*(by(i+1)+by(i))*(by(i+1)-by(i))/2)
     &        /dxm(i)/((ro(i)+ro(i+1))/2))
      enddo

      do i=2,ix
        du=vxm(i)-vxm(i-1)
        roh(i)=roh(i)+dt*(-du/dx(i)*ro(i))

        du=(vxm(i)-vxm(i-1)+vxmh(i)-vxmh(i-1))/2
        teh(i)=teh(i)+dt*(-du/dx(i)*(gm-1)*gm*pr(i)/ro(i))
      enddo

      call cipdxsrc(rodxh,ro,roh,vx,dt,dx,ix)
      call cipdxsrc(tedxh,te,teh,vx,dt,dx,ix)
      call cipdxsrc(vydxh,vy,vyh,vx,dt,dx,ix)
      call cipdxsrc(vxdxmh,vxm,vxmh,vxm,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call cipadv(roh,rodxh,vx,0,dt,dxm,ix)
      call cipadv(teh,tedxh,vx,0,dt,dxm,ix)
      call cipadv(vyh,vydxh,vx,0,dt,dxm,ix)
      call cipadv(vxmh,vxdxmh,vxm,1,dt,dx,ix)

c----------------------------------------------------------------------|
c--- MOC/CT Method
c----------------------------------------------------------------------|
c--- ### important ### Use (n+1)-step values for vx !

      call moc(vyc,byc,ro,vxmh,bxm,vyh,by,dt,dxm,ix)
      call ctranspt(byh,vxmh,vyc,bxm,byc,dt,dx,ix)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do i=1,ix
        ro(i)=roh(i)
        te(i)=teh(i)
        vxm(i)=vxmh(i)
        vy(i)=vyh(i)
        rodx(i)=rodxh(i)
        tedx(i)=tedxh(i)
        vxdxm(i)=vxdxmh(i)
        vydx(i)=vydxh(i)
        by(i)=byh(i)
      enddo

c     calculate pressure
      do i=1,ix
         pr(i) = ro(i)*te(i)/gm
      enddo

      do i=2,ix-1
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo

      return
      end
