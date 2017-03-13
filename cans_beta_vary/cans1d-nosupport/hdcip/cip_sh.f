c======================================================================|
      subroutine cip_sh(de,ei,rxm
     &   ,dedx,eidx,rxdxm,dt,cvis,gm,dx,dxm,ix)
c======================================================================|
c
c NAME  cip_sh
c
c PURPOSE
c    solve eqs. by CIP method
c        * hydrodynamics
c        * special relativity
c
c INPUTS & OUTPUTS
c    ro(ix): [double] density
c    ei(ix): [double] internal energy density
c    vxm(ix): [double] velocity 
c    rodx(ix): [double] density gradient
c    eidx(ix): [double] int-energy  gradient
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

      dimension  ro(ix),pr(ix),vx(ix),vxm(ix)
      dimension  de(ix),ei(ix),rxm(ix)
      dimension  gl(ix),glm(ix)
      dimension  deh(ix),eih(ix),rxmh(ix)
      dimension  dedx(ix),eidx(ix),rxdxm(ix)
      dimension  dedxh(ix),eidxh(ix),rxdxmh(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do i=1,ix
        deh(i)=de(i)
        eih(i)=ei(i)
        rxmh(i)=rxm(i)
        dedxh(i)=dedx(i)
        eidxh(i)=eidx(i)
        rxdxmh(i)=rxdxm(i)
      enddo

      do i=1,ix-1
        dem=(de(i)+de(i+1))/2
        eim=(ei(i)+ei(i+1))/2
        rtm = sqrt( rxm(i)**2 + ( dem+gm*eim )**2 )
        vxm(i) = rxm(i)/rtm
        glm(i) = 1.d0/sqrt(1.-vxm(i)**2)
      enddo
      do i=2,ix
        rx=(rxm(i-1)+rxm(i))/2
        rt = sqrt( rx**2 + ( de(i)+gm*ei(i) )**2 )
        vx(i) = rx/rt
        gl(i) = 1./sqrt(1.-vx(i)**2)
        pr(i) = ei(i)/gl(i)*(gm-1.d0)
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|

c     cvis=0.75
c     cvis=3.00
      do i=2,ix
        du=vxm(i)-vxm(i-1)
        if (du.lt.0.0) then
          cs=sqrt(ei(i)/de(i))
          pr(i)=pr(i)+cvis*(-de(i)*cs*du+(gm+1)/2*de(i)*du**2)
        endif
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### vx must be advanced first.

      do i=1,ix-1
        rxmh(i)=rxmh(i)+dt*(-(pr(i+1)-pr(i))/dxm(i))
     &                 +dt*(-rxm(i)*(vx(i+1)-vx(i))/dxm(i))
      enddo

      do i=2,ix
        du=vxm(i)-vxm(i-1)
        deh(i)=deh(i)+dt*(-du/dx(i)*de(i))

        du=vxm(i)-vxm(i-1)
        dg=glm(i)*vxm(i)-glm(i-1)*vxm(i-1)
        eih(i)=eih(i)+dt*(-du/dx(i)*ei(i))
     &               +dt*(-dg/dx(i)*pr(i))
      enddo

      call cipdxsrc(dedxh,de,deh,vx,dt,dx,ix)
      call cipdxsrc(eidxh,ei,eih,vx,dt,dx,ix)
      call cipdxsrc(rxdxmh,rxm,rxmh,vxm,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call cipadv(deh,dedxh,vx,0,dt,dxm,ix)
      call cipadv(eih,eidxh,vx,0,dt,dxm,ix)
      call cipadv(rxmh,rxdxmh,vxm,1,dt,dx,ix)

c----------------------------------------------------------------------|
c--- time derivative term ---
c----------------------------------------------------------------------|
      do i=2,ix
        rxh=(rxmh(i-1)+rxmh(i))/2
        rth = sqrt( rxh**2 + ( deh(i)+gm*eih(i) )**2 )
        vxh = rxh/rth
        glh = 1./sqrt(1.-vxh**2)
        h=(gm-1.)*(glh-gl(i))/((glh+gl(i))/2)
        eih(i)=eih(i)*exp(-h)
      enddo

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|

      do i=1,ix
        de(i)=deh(i)
        ei(i)=eih(i)
        rxm(i)=rxmh(i)
        dedx(i)=dedxh(i)
        eidx(i)=eidxh(i)
        rxdxm(i)=rxdxmh(i)
      enddo


      return
      end
