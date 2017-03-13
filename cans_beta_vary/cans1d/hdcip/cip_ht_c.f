c======================================================================|
      subroutine cip_ht_c(ro,vx,vxm,rodx,vxdxm
     &   ,dt,cvis,cs2,sc,scm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_ht_c
c
c PURPOSE
c    solve eqs. by CIP  method with effects of
c        * isothermal hydrodynamics
c        * non-uniform cross section
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    vxm(ix): [double] velocity 
c    rodx(ix): [double] density gradient
c    vxdxm(ix): [double] velocity gradient
c
c OUTPUTS
c    vx(ix): [double] velocity 
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    sc(ix), scm(ix) : [double] cross section
c    dsc(ix), dscm(ix) : [double] cross section gradient
c    dx(ix), dxm(ix) : [double] grid spacing
c    cs2: [double] square of sound speed
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix),vxm(ix)
      dimension roh(ix),vxmh(ix)
      dimension rodx(ix),vxdxm(ix)
      dimension rodxh(ix),vxdxmh(ix)
      dimension sc(ix),scm(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix
        roh(i)=ro(i)
        vxmh(i)=vxm(i)
        rodxh(i)=rodx(i)
        vxdxmh(i)=vxdxm(i)
      enddo

c--- velocity ---
      do i=2,ix
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo

      do i=1,ix
         pr(i) = cs2*ro(i)
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|

c     cvis=0.75
      do i=2,ix
        du=vxm(i)-vxm(i-1)
        if (du.lt.0.0) then
          cs=sqrt(cs2)
          pr(i)=pr(i)+cvis*(-ro(i)*cs*du+ro(i)*du**2)
        endif
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### vx must be advanced first.

      do i=1,ix-1
        vxmh(i)=vxmh(i)+dt*(-(pr(i+1)-pr(i))/dxm(i)/((ro(i)+ro(i+1))/2))
      enddo

      do i=2,ix
        du=vxm(i)*scm(i)-vxm(i-1)*scm(i-1)
        roh(i)=roh(i)+dt*(-du/dx(i)*ro(i)/sc(i))
      enddo

      call cipdxsrc(rodxh,ro,roh,vx,dt,dx,ix)
      call cipdxsrc(vxdxmh,vxm,vxmh,vxm,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call cipadv(roh,rodxh,vx,0,dt,dxm,ix)
      call cipadv(vxmh,vxdxmh,vxm,1,dt,dx,ix)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|

      do i=1,ix
        ro(i)=roh(i)
        vxm(i)=vxmh(i)
        rodx(i)=rodxh(i)
        vxdxm(i)=vxdxmh(i)
      enddo

c     calculate pressure

      do i=2,ix-1
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo


      return
      end
