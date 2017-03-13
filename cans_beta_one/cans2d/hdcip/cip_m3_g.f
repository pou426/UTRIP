c======================================================================|
      subroutine cip_m3_g(ro,pr,vx,vxm,vy,vym,vz,te
     &             ,rodx,rody,tedx,tedy
     &             ,vxdxm,vxdym,vydxm,vydym,vzdx,vzdy
     &             ,bx,bxm,by,bym,bz,az
     &             ,dt,cvis,gm,gxm,gym,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  cip_m
c
c PURPOSE
c    solve eqs. by CIP method with effects of
c        * MHD
c
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    rodx(ix,jx): [double] density gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    vx(ix,jx), vxm(ix,jx) : [double] velocity along the x-cordinate
c    vy(ix,jx), vym(ix,jx) : [double] velocity along the y-cordinate
c    dx(ix), dxm(ix) : [double] grid spacing
c    dy(jx), dym(jx) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension ro(ix,jx),te(ix,jx),vxm(ix,jx),vym(ix,jx)
      dimension vz(ix,jx)
      dimension rodx(ix,jx),rody(ix,jx),tedx(ix,jx),tedy(ix,jx)
      dimension vxdxm(ix,jx),vxdym(ix,jx),vydxm(ix,jx),vydym(ix,jx)
      dimension vzdx(ix,jx),vzdy(ix,jx)
      dimension roh(ix,jx),teh(ix,jx),vxmh(ix,jx),vymh(ix,jx)
      dimension vzh(ix,jx)
      dimension rodxh(ix,jx),rodyh(ix,jx),tedxh(ix,jx),tedyh(ix,jx)
      dimension vzdxh(ix,jx),vzdyh(ix,jx)
      dimension vxdxmh(ix,jx),vxdymh(ix,jx),vydxmh(ix,jx),vydymh(ix,jx)
      dimension pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension vxny(ix,jx),vynx(ix,jx)

      dimension bx(ix,jx),bxm(ix,jx),bxmh(ix,jx)
      dimension by(ix,jx),bym(ix,jx),bymh(ix,jx)
      dimension bz(ix,jx),bzh(ix,jx)
      dimension vxc(ix,jx),vyc(ix,jx)
      dimension bxc(ix,jx),byc(ix,jx)
      dimension bznx(ix,jx),bzny(ix,jx)
      dimension vznx(ix,jx),vzny(ix,jx)
      dimension az(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4

c---  calculate specific thermal energy from pressure

      do j=1,jx
      do i=1,ix
        roh(i,j)=ro(i,j)
        rodxh(i,j)=rodx(i,j)
        rodyh(i,j)=rody(i,j)
        teh(i,j)=te(i,j)
        tedxh(i,j)=tedx(i,j)
        tedyh(i,j)=tedy(i,j)
        vxmh(i,j)=vxm(i,j)
        vxdxmh(i,j)=vxdxm(i,j)
        vxdymh(i,j)=vxdym(i,j)
        vymh(i,j)=vym(i,j)
        vydxmh(i,j)=vydxm(i,j)
        vydymh(i,j)=vydym(i,j)
        bxmh(i,j) = bxm(i,j)
        bymh(i,j) = bym(i,j)
        vzh(i,j)=vz(i,j)
        vzdxh(i,j)=vzdx(i,j)
        vzdyh(i,j)=vzdy(i,j)
        bzh(i,j) = bz(i,j)
      enddo
      enddo

      do j=1,jx
      do i=1,ix
        pr(i,j) = ro(i,j)*te(i,j)/gm
      enddo
      enddo

      do j=1,jx
      do i=2,ix
        vx(i,j) = (vxm(i-1,j)+vxm(i,j))/2
        bx(i,j) = (bxm(i-1,j)+bxm(i,j))/2
      enddo
      enddo
      do j=2,jx
      do i=1,ix
        vy(i,j) = (vym(i,j-1)+vym(i,j))/2
        by(i,j) = (bym(i,j-1)+bym(i,j))/2
      enddo
      enddo
      do j=1,jx-1
      do i=2,ix
        vxny(i,j) = (vxm(i-1,j)+vxm(i,j)+vxm(i-1,j+1)+vxm(i,j+1))/4
      enddo
      enddo
      do j=2,jx
      do i=1,ix-1
        vynx(i,j) = (vym(i,j-1)+vym(i,j)+vym(i+1,j-1)+vym(i+1,j))/4
      enddo
      enddo
      do j=2,jx-1
      do i=2,ix-1
         az(i,j)=az(i,j)-dt*(-vx(i,j)*by(i,j)+vy(i,j)*bx(i,j))
      enddo
      enddo
c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|
c     cvis=0.75
      do j=2,jx
      do i=2,ix
        dvx=vxm(i,j)-vxm(i-1,j)
        dvy=vym(i,j)-vym(i,j-1)
        divv=dvx/dx(i)+dvy/dy(j)
        if (divv.lt.0.0) then
          du=divv*min(dx(i),dy(j))
          cs=sqrt(abs(te(i,j)))
          pr(i,j)=pr(i,j)
     &            +cvis*(-ro(i,j)*cs*du+(gm+1)/2*ro(i,j)*du**2)
        endif
      enddo
      enddo
c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### Calculate magnetic stress term Using MOC method

      call moclag(bxc,byc,ro,vxm,bxm,vym,bym,dt,dxm,dym,ix,jx)
      call moclagz(bznx,bzny,ro,bxm,bym,vz,bz,dt,dxm,dym,ix,jx)
      do j=1,jx-1
      do i=2,ix-1
        bxny=(bxm(i,j)+bxm(i-1,j)+bxm(i,j+1)+bxm(i-1,j+1))/4
        rony=(ro(i,j)+ro(i,j+1))/2
        vymh(i,j)=vymh(i,j)
     &    +1./(4.*pi*rony)
     &      *bxny*(byc(i,j)-byc(i-1,j))/dx(i)*dt
      enddo
      enddo
      do j=2,jx-1
      do i=1,ix-1
        bynx=(bym(i,j)+bym(i+1,j)+bym(i,j-1)+bym(i+1,j-1))/4
        ronx=(ro(i,j)+ro(i+1,j))/2
        vxmh(i,j)=vxmh(i,j)
     &    +1./(4.*pi*ronx)
     &      *bynx*(bxc(i,j)-bxc(i,j-1))/dy(j)*dt
      enddo
      enddo
      do j=2,jx-1
      do i=2,ix-1
        vzh(i,j)=vzh(i,j)
     &    +1./(4.*pi*ro(i,j))
     &      *(bx(i,j)*(bznx(i,j)-bznx(i-1,j))/dx(i)*dt
     &       +by(i,j)*(bzny(i,j)-bzny(i,j-1))/dy(j)*dt)
      enddo
      enddo

c--- ### important ### velocity must be advanced first.

      do j=1,jx
      do i=1,ix-1
        vxmh(i,j)=vxmh(i,j)-2*dt/dxm(i)/(ro(i,j)+ro(i+1,j))
     &    *(pr(i+1,j)-pr(i,j)
     &     +pi4i*(by(i+1,j)**2+bz(i+1,j)**2-by(i,j)**2-bz(i,j)**2)/2)
     &       +dt*gxm(i,j)
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
        vymh(i,j)=vymh(i,j)-2*dt/dym(j)/(ro(i,j)+ro(i,j+1))
     &    *(pr(i,j+1)-pr(i,j)
     &     +pi4i*(bx(i,j+1)**2+bz(i,j+1)**2-bx(i,j)**2-bz(i,j)**2)/2)
     &       +dt*gym(i,j)
      enddo
      enddo

      do j=2,jx
      do i=2,ix
        dvx=vxm(i,j)-vxm(i-1,j)
        dvy=vym(i,j)-vym(i,j-1)
        divv=dvx/dx(i)+dvy/dy(j)
        roh(i,j)=roh(i,j)+dt*(-divv*ro(i,j))

        dvx=(vxm(i,j)-vxm(i-1,j)+vxmh(i,j)-vxmh(i-1,j))/2
        dvy=(vym(i,j)-vym(i,j-1)+vymh(i,j)-vymh(i,j-1))/2
        divv=dvx/dx(i)+dvy/dy(j)
        teh(i,j)=teh(i,j)+dt*(-divv*(gm-1)*gm*pr(i,j)/ro(i,j))
      enddo
      enddo


      call cipdxsrc(rodxh,rodyh,ro,roh,vx,vy,dt,dx,dy,ix,jx)
      call cipdxsrc(tedxh,tedyh,te,teh,vx,vy,dt,dx,dy,ix,jx)
      call cipdxsrc(vxdxmh,vxdymh,vxm,vxmh,vxm,vynx,dt,dx,dy,ix,jx)
      call cipdxsrc(vydxmh,vydymh,vym,vymh,vxny,vym,dt,dx,dy,ix,jx)
      call cipdxsrc(vzdxh,vzdyh,vz,vzh,vx,vy,dt,dx,dy,ix,jx)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      call cipadv(roh,rodxh,rodyh,vx,vy,0,0,dt,dxm,dym,ix,jx)
      call cipadv(teh,tedxh,tedyh,vx,vy,0,0,dt,dxm,dym,ix,jx)
      call cipadv(vzh,vzdxh,vzdyh,vx,vy,0,0,dt,dxm,dym,ix,jx)
      call cipadv(vxmh,vxdxmh,vxdymh,vxm,vynx,1,0,dt,dx,dym,ix,jx)
      call cipadv(vymh,vydxmh,vydymh,vxny,vym,0,1,dt,dxm,dy,ix,jx)

c----------------------------------------------------------------------|
c--- MOC/CT Method
c----------------------------------------------------------------------|
c--- ### important ### Use (n+1)-step values for vx & vy !
      call moc(vxc,vyc,bxc,byc,ro,vxmh,bxm,vymh,bym,dt,dxm,dym,ix,jx)
      call ctranspt(bxmh,bymh,vxc,vyc,bxc,byc,dt,dx,dy,ix,jx)

      call mocz(vznx,vzny,bznx,bzny,ro,vxmh,bxm,vymh,bym,vz,bz
     &       ,dt,dxm,dym,ix,jx)
      call ctransptz(bzh,vznx,vzny,bznx,bzny,vxmh,vymh,bxm,bym
     &       ,dt,dx,dy,ix,jx)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
        ro(i,j)=roh(i,j)
        rodx(i,j)=rodxh(i,j)
        rody(i,j)=rodyh(i,j)
        te(i,j)=teh(i,j)
        tedx(i,j)=tedxh(i,j)
        tedy(i,j)=tedyh(i,j)
        vxm(i,j)=vxmh(i,j)
        vxdxm(i,j)=vxdxmh(i,j)
        vxdym(i,j)=vxdymh(i,j)
        vym(i,j)=vymh(i,j)
        vydxm(i,j)=vydxmh(i,j)
        vydym(i,j)=vydymh(i,j)
        vz(i,j)=vzh(i,j)
        vzdx(i,j)=vzdxh(i,j)
        vzdy(i,j)=vzdyh(i,j)
        bxm(i,j)=bxmh(i,j)
        bym(i,j)=bymh(i,j)
        bz(i,j)=bzh(i,j)
      enddo
      enddo

      do j=1,jx
      do i=1,ix
        pr(i,j) = ro(i,j)*te(i,j)/gm
      enddo
      enddo

      do j=1,jx
      do i=2,ix
        vx(i,j) = (vxm(i-1,j)+vxm(i,j))/2
        bx(i,j) = (bxm(i-1,j)+bxm(i,j))/2
      enddo
      enddo
      do j=2,jx
      do i=1,ix
        vy(i,j) = (vym(i,j-1)+vym(i,j))/2
        by(i,j) = (bym(i,j-1)+bym(i,j))/2
      enddo
      enddo

      return
      end
