c======================================================================|
      subroutine ccip_h(rdm,ro,pr,vx,te,vxm,rodx,tedx,vxdxm
     &     ,dt,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_h
c
c PURPOSE
c    solve eqs. by conservative CIP method with effects of
c        * hydrodynamics
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    te(ix): [double] temperature
c    vxm(ix): [double] velocity 
c    rodx(ix): [double] density gradient
c    tedx(ix): [double] temperature gradient
c    vxdxm(ix): [double] velocity gradient
c
c OUTPUTS
c    vx(ix): [double] velocity 
c    pr(ix): [double] pressure
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
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

      dimension  ro(ix),pr(ix),te(ix),vx(ix),vxm(ix)
      dimension  roh(ix),teh(ix),vxmh(ix)
      dimension  rodx(ix),tedx(ix),vxdxm(ix)
      dimension  rodxh(ix),tedxh(ix),vxdxmh(ix)
      dimension  rdm(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix
        roh(i)=ro(i)
        teh(i)=te(i)
        vxmh(i)=vxm(i)
        rodxh(i)=rodx(i)
        tedxh(i)=tedx(i)
        vxdxmh(i)=vxdxm(i)
      enddo

c     calculate pressure
      do i=1,ix
         pr(i) = ro(i)*te(i)/gm
      enddo
c--- velocity ---
      do i=2,ix
        vx(i)=(vxm(i)+vxm(i-1))/2
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|

      cvis=0.75
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
c--- ### important ### vx must be advanced first.

      do i=1,ix-1
        vxmh(i)=vxmh(i)+dt*(-(pr(i+1)-pr(i))/dxm(i)/((ro(i)+ro(i+1))/2))
      enddo

      do i=2,ix
        du=vxm(i)-vxm(i-1)
        roh(i)=roh(i)+dt*(-du/dx(i)*ro(i))

        du=(vxm(i)-vxm(i-1)+vxmh(i)-vxmh(i-1))/2
        teh(i)=teh(i)+dt*(-du/dx(i)*(gm-1)*gm*pr(i)/ro(i))
      enddo

      call cipdxsrc(rodxh,ro,roh,vx,dt,dx,ix)
      call cipdxsrc(tedxh,te,teh,vx,dt,dx,ix)
      call cipdxsrc(vxdxmh,vxm,vxmh,vxm,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call ccipadvrd(rdm,roh,dt,vx,dxm,ix)
      call cipadv(teh,tedxh,vx,0,dt,dxm,ix)
      call cipadv(vxmh,vxdxmh,vxm,1,dt,dx,ix)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|

      do i=1,ix
        ro(i)=roh(i)
        te(i)=teh(i)
        vxm(i)=vxmh(i)
        rodx(i)=rodxh(i)
        tedx(i)=tedxh(i)
        vxdxm(i)=vxdxmh(i)
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
