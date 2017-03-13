c======================================================================|
      subroutine cip_mt_bg(ro,vx,vy,by,vxm,rodx,vxdxm,vyrdx
     &    ,bx,bxm,dt,cvis,cs2,gxm,sc,scm,rr,rrm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_mt_bg
c
c PURPOSE
c    solve eqs. by CIP-MOCCT  method with effects of
c        * isothermal MHD
c        * axial symmetry
c        * non-uniform poloidal magnetic field
c        * gravity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vxm(ix): [double] velocity 
c    vy(ix): [double] velocity 
c    vz(ix): [double] velocity 
c    by(ix): [double] magnetic field
c    bz(ix): [double] magnetic field
c    rodx(ix): [double] density gradient
c    vxdxm(ix): [double] velocity gradient
c    vyrdx(ix): [double] velocity gradient
c    vzdx(ix): [double] velocity gradient
c
c OUTPUTS
c    vx(ix): [double] velocity 
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    bx(ix), bxm(ix) : [double] magnetic field
c    gx(ix), gxm(ix) : [double] gravity
c    gy(ix), gym(ix) : [double] gravity
c    gz(ix), gzm(ix) : [double] gravity
c    sc(ix), scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    rr(ix), rrm(ix) : [double] distance from rotation axis
c    drr(ix), drrm(ix) : [double] gradient of distance from rotation axis
c    dx(ix), dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix),vyr(ix)
      dimension vxm(ix)
      dimension bx(ix),bxm(ix)
      dimension roh(ix),vxmh(ix),vyrh(ix)
      dimension rodx(ix),vxdxm(ix),vyrdx(ix)
      dimension rodxh(ix),vxdxmh(ix),vyrdxh(ix)
      dimension vym(ix),bym(ix)
      dimension sc(ix),scm(ix)
      dimension gxm(ix)
      dimension rr(ix),rrm(ix)
      dimension byscm(ix),bxscm(ix),bysch(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c---  calculate specific thermal energy from pressure

      do i=1,ix
        roh(i) = ro(i)
        vxmh(i) = vxm(i)
        vyr(i)  = vy(i)*rr(i)
        vyrh(i) = vyr(i)
        rodxh(i)=rodx(i)
        vxdxmh(i)=vxdxm(i)
        vyrdxh(i)=vyrdx(i)
        bysch(i)=by(i)*sc(i)/rr(i)
      enddo

c     calculate pressure
      do i=1,ix
         pr(i) = ro(i)*cs2
      enddo
c--- velocity on the grid points ---
      do i=2,ix
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|

c     cvis=0.75
      cs=sqrt(cs2)
      do i=2,ix
        du=vxm(i)-vxm(i-1)
        if (du.lt.0.0) then
          pr(i)=pr(i)+cvis*(-ro(i)*cs*du+ro(i)*du**2)
        endif
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### Calculate magnetc stress term Using MOC method

      call moclag(bym,ro,bxm,vy,by,dt,dxm,ix)
      do i=2,ix-1
        vyrh(i)=vyrh(i)
     &    +1./(4.*pi*ro(i))*bx(i)
     &    *(rrm(i)*bym(i)-rrm(i-1)*bym(i-1))/dx(i)*dt
      enddo

c--- ### important ### vx must be advanced first.
      do i=1,ix-1
        vxmh(i)=vxmh(i)+dt*(
     &    -( pr(i+1)-pr(i) )/dxm(i) /((ro(i)+ro(i+1))/2)
     &    -( rr(i+1)*by(i+1)-rr(i)*by(i) )/dxm(i)
     &       *pi4i*(by(i+1)+by(i))/2/rrm(i) /((ro(i)+ro(i+1))/2) 
     &    +( rr(i+1)-rr(i) )/dxm(i)
     &       *((vy(i+1)+vy(i))/2)**2/rrm(i)
     &     )
     &          +dt*gxm(i)
      enddo

      do i=2,ix
        du=vxm(i)*scm(i)-vxm(i-1)*scm(i-1)
        roh(i)=roh(i)+dt*(-du/dx(i)*ro(i)/sc(i))
      enddo

      call cipdxsrc(rodxh,ro,roh,vx,dt,dx,ix)
      call cipdxsrc(vyrdxh,vyr,vyrh,vx,dt,dx,ix)
      call cipdxsrc(vxdxmh,vxm,vxmh,vxm,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call cipadv(roh,rodxh,vx,0,dt,dxm,ix)
      call cipadv(vyrh,vyrdxh,vx,0,dt,dxm,ix)
      call cipadv(vxmh,vxdxmh,vxm,1,dt,dx,ix)

c----------------------------------------------------------------------|
c--- MOC/CT Method
c----------------------------------------------------------------------|
c--- ### important ### Use (n+1)-step values for vx !

      call moc(vym,bym,ro,vxmh,bxm,vy,by,dt,dxm,ix)
      do i=1,ix
        bxscm(i)=bxm(i)*scm(i)/rrm(i)
        byscm(i)=bym(i)*scm(i)/rrm(i)
      enddo
      call ctranspt(bysch,vxmh,vym,bxscm,byscm,dt,dx,ix)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do i=1,ix
        ro(i)=roh(i)
        vxm(i)=vxmh(i)
        vy(i)=vyrh(i)/rr(i)
        rodx(i)=rodxh(i)
        vxdxm(i)=vxdxmh(i)
        vyrdx(i)=vyrdxh(i)
        by(i)=bysch(i)/sc(i)*rr(i)
      enddo

c     calculate pressure
      do i=2,ix-1
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo


      return
      end
