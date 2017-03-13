c======================================================================|
      subroutine cip_sm(de,ei,rxm,ry,by,dedx,eidx,rxdxm,rydx
     &    ,bx,bxm,ro,pr,dt,cvis,gm,dx,dxm,ix)
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

      dimension  by(ix),ry(ix),vy(ix)
      dimension  bx(ix),bxm(ix)
      dimension  byh(ix),ryh(ix),vyh(ix)
      dimension  rydx(ix),rydxh(ix)
      dimension  vyc(ix),byc(ix)
      dimension  pm2(ix)
      dimension  bbt(ix),bbx(ix),bby(ix)
      dimension  bbtm(ix),bbxm(ix),bbym(ix)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1.d0/pi/4.d0
      pi4ir=1.d0/sqrt(pi)/2.d0

      do i=1,ix
        deh(i)=de(i)
        eih(i)=ei(i)
        rxmh(i)=rxm(i)
        dedxh(i)=dedx(i)
        eidxh(i)=eidx(i)
        rxdxmh(i)=rxdxm(i)
        ryh(i) = ry(i)
        rydxh(i)=rydx(i)
        byh(i) = by(i)
      enddo

      do i=1,ix-1
        dem=(de(i)+de(i+1))/2
        eim=(ei(i)+ei(i+1))/2
        rym=(ry(i)+ry(i+1))/2
        bym=(by(i)+by(i+1))/2
        rtm = sqrt( rxm(i)**2 + rym**2 + ( dem+gm*eim )**2 )
        vxm(i) = rxm(i)/rtm
        vym    = rym   /rtm
        glm(i) = 1.d0/sqrt(1.-vxm(i)**2-vym**2)
        bbtm(i)=glm(i)*pi4ir*(vxm(i)*bxm(i)+vym*bym)
        bbxm(i)=bxm(i)/glm(i)*pi4ir+vxm(i)*bbtm(i)
        bbym(i)=bym   /glm(i)*pi4ir+vym   *bbtm(i)
      enddo
      do i=1,ix
        gl(i)=1.d0
        vx(i)=0.d0
        vy(i)=0.d0
      enddo
      do i=2,ix
        rx=(rxm(i-1)+rxm(i))/2
        rt = sqrt( rx**2 + ry(i)**2 + ( de(i)+gm*ei(i) )**2 )
        vx(i) = rx   /rt
        vy(i) = ry(i)/rt
        vyh(i) = vy(i)
        gl(i) = 1./sqrt(1.-vx(i)**2-vy(i)**2)
        bbt(i)=gl(i)*pi4ir*(vx(i)*bx(i)+vy(i)*by(i))
        bbx(i)=bx(i)/gl(i)*pi4ir+vx(i)*bbt(i)
        bby(i)=by(i)/gl(i)*pi4ir+vy(i)*bbt(i)
      enddo

      do i=1,ix-1
        pm2(i)=pi4i*(bx(i)**2+by(i)**2)/gl(i)**2
     &        +pi4i*(vx(i)*bx(i)+vy(i)*by(i))**2
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
c--- ### important ### Calculate magnetc stress term Using MOC method

      call moclag_sm(byc,ro,bxm,vy,by,pr,glm,pm2,gm,dt,dxm,ix)
      do i=2,ix-1
        bbt(i)=gl(i)*pi4ir*(vx(i)*bx(i)+vy(i)*byc(i))
        bbx(i)=bx(i)/gl(i)*pi4ir+vx(i)*bbt(i)
        bby(i)=byc(i)/gl(i)*pi4ir+vy(i)*bbt(i)
        ryh(i)=ryh(i)+bbx(i)*(bby(i)-bby(i-1))/dx(i)*dt
     &    +dt*(-ry(i)*(vxm(i)-vxm(i-1))/dx(i))
      enddo

c--- ### important ### vx must be advanced first.

      do i=1,ix-1
        rxmh(i)=rxmh(i)
     &    +dt*(-(pr(i+1)-pr(i)+pm2(i+1)/2-pm2(i)/2
     &         -bbxm(i+1)**2+bbxm(i)**2)/dxm(i))
     &    +dt*(-rxm(i)*(vx(i+1)-vx(i))/dxm(i))
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
      call cipdxsrc(rydxh,ry,ryh,vx,dt,dx,ix)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|
      call cipadv(deh,dedxh,vx,0,dt,dxm,ix)
      call cipadv(eih,eidxh,vx,0,dt,dxm,ix)
      call cipadv(rxmh,rxdxmh,vxm,1,dt,dx,ix)
      call cipadv(ryh,rydxh,vx,0,dt,dxm,ix)

c----------------------------------------------------------------------|
c--- MOC/CT Method
c----------------------------------------------------------------------|
c--- ### important ### Use (n+1)-step values for vx !

      call moc_sm(vyc,byc,ro,vxm,bxm,vyh,by,pr,glm,pm2,gm,dt,dxm,ix)
      call ctranspt(byh,vxm,vyc,bxm,byc,dt,dx,ix)

c----------------------------------------------------------------------|
c--- time derivative term ---
c----------------------------------------------------------------------|
      do i=1,ix-1
        dem=(deh(i)+deh(i+1))/2
        eim=(eih(i)+eih(i+1))/2
        rym=(ryh(i)+ryh(i+1))/2
        bym=(byh(i)+byh(i+1))/2
        rtm = sqrt( rxmh(i)**2 + rym**2 + ( dem+gm*eim )**2 )
        vxm(i) = rxmh(i)/rtm
        vym    = rym   /rtm
        glm(i) = 1.d0/sqrt(1.-vxm(i)**2-vym**2)
        bbtmh   =glm(i)*pi4ir*(vxm(i)*bxm(i)+vym*bym)
        bbxmh   =bxm(i)/glm(i)*pi4ir+vxm(i)*bbtmh
        bbymh   =bym   /glm(i)*pi4ir+vym   *bbtmh
        rxmh(i)=rxmh(i)+(bbxmh   *bbtmh   -bbxm(i)*bbtm(i))
      enddo
      do i=2,ix
        rxh=(rxmh(i-1)+rxmh(i))/2
        rth = sqrt( rxh**2 + ryh(i)**2 + ( deh(i)+gm*eih(i) )**2 )
        vxh = rxh/rth
        vyh(i) = ryh(i)/rth
        glh = 1./sqrt(1.-vxh**2-vyh(i)**2)
        h=(gm-1.)*(glh-gl(i))/((glh+gl(i))/2)
        eih(i)=eih(i)*exp(-h)

        bbth  =glh  *pi4ir*(vxh*bx(i)+vyh(i)*byh(i))
        bbyh  =byh(i)/glh*pi4ir+vyh(i)*bbth
        ryh(i)=ryh(i)+(bbyh*bbth-bby(i)*bbt(i))
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
        ry(i) = ryh(i)
        rydx(i)=rydxh(i)
        by(i) = byh(i)
      enddo


      return
      end
