c======================================================================|
      subroutine cip_m(ro,pr,vx,vxm,vy,vym,vz,vzm,te
     &                   ,rodx,rody,rodz,tedx,tedy,tedz
     &                   ,vxdxm,vxdym,vxdzm,vydxm,vydym,vydzm
     &                   ,vzdxm,vzdym,vzdzm,bx,bxm,by,bym,bz,bzm
     &                   ,dt,cvis,gm,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  cip_m
c
c PURPOSE
c    solve eqs. by CIP method with effects of
c        * MHD
c
c INPUTS & OUTPUTS
c    ro(ix,jx,kx): [double] density
c    rodx(ix,jx,kx): [double] density gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    vx(ix,jx,kx),vxm(ix,jx,kx): [double] velocity along the x-cordinate
c    vy(ix,jx,kx),vym(ix,jx,kx): [double] velocity along the y-cordinate
c    vz(ix,jx,kx),vzm(ix,jx,kx): [double] velocity along the z-cordinate
c    dx(ix), dxm(ix) : [double] grid spacing
c    dy(jx), dym(jx) : [double] grid spacing
c    dz(kx), dzm(kx) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c    modified 2004-11-1 K. Takahashi
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension ro(ix,jx,kx),te(ix,jx,kx)
      dimension vxm(ix,jx,kx),vym(ix,jx,kx),vzm(ix,jx,kx)
      dimension rodx(ix,jx,kx),rody(ix,jx,kx),rodz(ix,jx,kx)
      dimension tedx(ix,jx,kx),tedy(ix,jx,kx),tedz(ix,jx,kx)
      dimension vxdxm(ix,jx,kx),vxdym(ix,jx,kx),vxdzm(ix,jx,kx)
      dimension vydxm(ix,jx,kx),vydym(ix,jx,kx),vydzm(ix,jx,kx)
      dimension vzdxm(ix,jx,kx),vzdym(ix,jx,kx),vzdzm(ix,jx,kx)
      dimension roh(ix,jx,kx),teh(ix,jx,kx)
      dimension vxmh(ix,jx,kx),vymh(ix,jx,kx),vzmh(ix,jx,kx)
      dimension rodxh(ix,jx,kx),rodyh(ix,jx,kx),rodzh(ix,jx,kx)
      dimension tedxh(ix,jx,kx),tedyh(ix,jx,kx),tedzh(ix,jx,kx)
      dimension vxdxmh(ix,jx,kx),vxdymh(ix,jx,kx),vxdzmh(ix,jx,kx)
      dimension vydxmh(ix,jx,kx),vydymh(ix,jx,kx),vydzmh(ix,jx,kx)
      dimension vzdxmh(ix,jx,kx),vzdymh(ix,jx,kx),vzdzmh(ix,jx,kx)
      dimension pr(ix,jx,kx),vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension vxny(ix,jx,kx),vxnz(ix,jx,kx)
      dimension vynx(ix,jx,kx),vynz(ix,jx,kx)
      dimension vznx(ix,jx,kx),vzny(ix,jx,kx)

      dimension bx(ix,jx,kx),bxm(ix,jx,kx),bxmh(ix,jx,kx)
      dimension by(ix,jx,kx),bym(ix,jx,kx),bymh(ix,jx,kx)
      dimension bz(ix,jx,kx),bzm(ix,jx,kx),bzmh(ix,jx,kx)
      dimension vxcy(ix,jx,kx),vxcz(ix,jx,kx)
      dimension bxcy(ix,jx,kx),bxcz(ix,jx,kx)
      dimension vycx(ix,jx,kx),vycz(ix,jx,kx)
      dimension bycx(ix,jx,kx),bycz(ix,jx,kx)
      dimension vzcx(ix,jx,kx),vzcy(ix,jx,kx)
      dimension bzcx(ix,jx,kx),bzcy(ix,jx,kx)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c---  calculate specific thermal energy from pressure

      do k=1,kx
      do j=1,jx
      do i=1,ix
        roh(i,j,k)=ro(i,j,k)
        rodxh(i,j,k)=rodx(i,j,k)
        rodyh(i,j,k)=rody(i,j,k)
        rodzh(i,j,k)=rodz(i,j,k)
        teh(i,j,k)=te(i,j,k)
        tedxh(i,j,k)=tedx(i,j,k)
        tedyh(i,j,k)=tedy(i,j,k)
        tedzh(i,j,k)=tedz(i,j,k)
        vxmh(i,j,k)=vxm(i,j,k)
        vxdxmh(i,j,k)=vxdxm(i,j,k)
        vxdymh(i,j,k)=vxdym(i,j,k)
        vxdzmh(i,j,k)=vxdzm(i,j,k)
        vymh(i,j,k)=vym(i,j,k)
        vydxmh(i,j,k)=vydxm(i,j,k)
        vydymh(i,j,k)=vydym(i,j,k)
        vydzmh(i,j,k)=vydzm(i,j,k)
        vzmh(i,j,k)=vzm(i,j,k)
        vzdxmh(i,j,k)=vzdxm(i,j,k)
        vzdymh(i,j,k)=vzdym(i,j,k)
        vzdzmh(i,j,k)=vzdzm(i,j,k)
        bxmh(i,j,k) = bxm(i,j,k)
        bymh(i,j,k) = bym(i,j,k)
        bzmh(i,j,k) = bzm(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
        pr(i,j,k) = ro(i,j,k)*te(i,j,k)/gm
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=2,ix
        vx(i,j,k) = (vxm(i-1,j,k)+vxm(i,j,k))/2
        bx(i,j,k) = (bxm(i-1,j,k)+bxm(i,j,k))/2
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx
      do i=1,ix
        vy(i,j,k) = (vym(i,j-1,k)+vym(i,j,k))/2
        by(i,j,k) = (bym(i,j-1,k)+bym(i,j,k))/2
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx
      do i=1,ix
        vz(i,j,k) = (vzm(i,j,k-1)+vzm(i,j,k))/2
        bz(i,j,k) = (bzm(i,j,k-1)+bzm(i,j,k))/2
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx-1
      do i=2,ix
        vxny(i,j,k) = 0.25*(vxm(i-1,j,k)+vxm(i-1,j+1,k)
     &                     +vxm(i  ,j,k)+vxm(i  ,j+1,k))
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=1,jx
      do i=2,ix
        vxnz(i,j,k) = 0.25*(vxm(i-1,j,k)+vxm(i-1,j,k+1)
     &                     +vxm(i  ,j,k)+vxm(i  ,j,k+1))
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx
      do i=1,ix-1
        vynx(i,j,k) = 0.25*(vym(i,j-1,k)+vym(i+1,j-1,k)
     &                     +vym(i,j  ,k)+vym(i+1,j  ,k))
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=2,jx
      do i=1,ix
        vynz(i,j,k) = 0.25*(vym(i,j-1,k)+vym(i,j-1,k+1)
     &                     +vym(i,j  ,k)+vym(i,j  ,k+1))
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx
      do i=1,ix-1
        vznx(i,j,k) = 0.25*(vzm(i,j,k-1)+vzm(i+1,j,k-1)
     &                     +vzm(i,j,k  )+vzm(i+1,j,k  ))
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx-1
      do i=1,ix
        vzny(i,j,k) = 0.25*(vzm(i,j,k-1)+vzm(i,j+1,k-1)
     &                     +vzm(i,j,k  )+vzm(i,j+1,k  ))
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c--- artifical additional pressure
c----------------------------------------------------------------------|
c     cvis=0.6
      dvis=cvis
      do k=2,kx
      do j=2,jx
      do i=2,ix
        dvx=vxm(i,j,k)-vxm(i-1,j,k)
        dvy=vym(i,j,k)-vym(i,j-1,k)
        dvz=vzm(i,j,k)-vzm(i,j,k-1)
        divv=dvx/dx(i)+dvy/dy(j)+dvz/dz(k)
        if (divv.lt.0.0) then
          du=divv*min(dx(i),dy(j),dz(k))
          cs=sqrt(abs(te(i,j,k)))
          pr(i,j,k)=pr(i,j,k)
     &             +cvis*(-ro(i,j,k)*cs*du)
     &             +dvis*(gm+1)/2*ro(i,j,k)*du**2
c        else
c          du=divv*min(dx(i),dy(j),dz(k))
c          cs=sqrt(abs(te(i,j,k)))
c          pr(i,j,k)=pr(i,j,k)
c     &             +cvis*(-ro(i,j,k)*cs*du)
        endif
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
c--- ### important ### Calculate magnetic stress term Using MOC method

      call moclagx(bxcy,bxcz,ro,bym,bzm,vx,bx,dt,dxm,dym,dzm,ix,jx,kx)
      call moclagy(bycx,bycz,ro,bxm,bzm,vy,by,dt,dxm,dym,dzm,ix,jx,kx)
      call moclagz(bzcx,bzcy,ro,bxm,bym,vz,bz,dt,dxm,dym,dzm,ix,jx,kx)
      do k=2,kx-1
      do j=2,jx-1
      do i=1,ix-1
        bynx=(bym(i,j,k)+bym(i+1,j,k)+bym(i,j-1,k)+bym(i+1,j-1,k))/4
        bznx=(bzm(i,j,k)+bzm(i+1,j,k)+bzm(i,j,k-1)+bzm(i+1,j,k-1))/4
        ronx=(ro(i,j,k)+ro(i+1,j,k))/2
        vxmh(i,j,k)=vxmh(i,j,k)
     &    +1./(4.*pi*ronx)
     &      *( bynx*(bxcy(i,j,k)-bxcy(i,j-1,k))/dy(j)*dt
     &        +bznx*(bxcz(i,j,k)-bxcz(i,j,k-1))/dz(k)*dt)
      enddo
      enddo
      enddo
      do k=2,kx-1
      do j=1,jx-1
      do i=2,ix-1
        bxny=(bxm(i,j,k)+bxm(i-1,j,k)+bxm(i,j+1,k)+bxm(i-1,j+1,k))/4
        bzny=(bzm(i,j,k)+bzm(i,j,k-1)+bzm(i,j+1,k)+bzm(i,j+1,k-1))/4
        rony=(ro(i,j,k)+ro(i,j+1,k))/2
        vymh(i,j,k)=vymh(i,j,k)
     &    +1./(4.*pi*rony)
     &      *( bxny*(bycx(i,j,k)-bycx(i-1,j,k))/dx(i)*dt
     &        +bzny*(bycz(i,j,k)-bycz(i,j,k-1))/dz(k)*dt)
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=2,jx-1
      do i=2,ix-1
        bxnz=(bxm(i,j,k)+bxm(i-1,j,k)+bxm(i,j,k+1)+bxm(i-1,j,k+1))/4
        bynz=(bym(i,j,k)+bym(i,j-1,k)+bym(i,j,k+1)+bym(i,j-1,k+1))/4
        ronz=(ro(i,j,k)+ro(i,j,k+1))/2
        vzmh(i,j,k)=vzmh(i,j,k)
     &    +1./(4.*pi*ronz)
     &      *( bxnz*(bzcx(i,j,k)-bzcx(i-1,j,k))/dx(i)*dt
     &        +bynz*(bzcy(i,j,k)-bzcy(i,j-1,k))/dy(j)*dt)
      enddo
      enddo
      enddo

c--- ### important ### velocity must be advanced first.

      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        vxmh(i,j,k)=vxmh(i,j,k)-2*dt/dxm(i)/(ro(i,j,k)+ro(i+1,j,k))
     &             *( pr(i+1,j,k)-pr(i,j,k)
     &               +pi4i*(by(i+1,j,k)**2-by(i,j,k)**2)/2
     &               +pi4i*(bz(i+1,j,k)**2-bz(i,j,k)**2)/2)
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        vymh(i,j,k)=vymh(i,j,k)-2*dt/dym(j)/(ro(i,j,k)+ro(i,j+1,k))
     &             *( pr(i,j+1,k)-pr(i,j,k)
     &               +pi4i*(bx(i,j+1,k)**2-bx(i,j,k)**2)/2
     &               +pi4i*(bz(i,j+1,k)**2-bz(i,j,k)**2)/2)
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        vzmh(i,j,k)=vzmh(i,j,k)-2*dt/dzm(k)/(ro(i,j,k)+ro(i,j,k+1))
     &             *( pr(i,j,k+1)-pr(i,j,k)
     &               +pi4i*(bx(i,j,k+1)**2-bx(i,j,k)**2)/2
     &               +pi4i*(by(i,j,k+1)**2-by(i,j,k)**2)/2)
      enddo
      enddo
      enddo

      do k=2,kx
      do j=2,jx
      do i=2,ix
        dvx=vxm(i,j,k)-vxm(i-1,j,k)
        dvy=vym(i,j,k)-vym(i,j-1,k)
        dvz=vzm(i,j,k)-vzm(i,j,k-1)
        divv=dvx/dx(i)+dvy/dy(j)+dvz/dz(k)
        roh(i,j,k)=roh(i,j,k)+dt*(-divv*ro(i,j,k))
        dvx=(vxm(i,j,k)-vxm(i-1,j,k)+vxmh(i,j,k)-vxmh(i-1,j,k))/2
        dvy=(vym(i,j,k)-vym(i,j-1,k)+vymh(i,j,k)-vymh(i,j-1,k))/2
        dvz=(vzm(i,j,k)-vzm(i,j,k-1)+vzmh(i,j,k)-vzmh(i,j,k-1))/2
        divv=dvx/dx(i)+dvy/dy(j)+dvz/dz(k)
        teh(i,j,k)=teh(i,j,k)+dt*(-divv*(gm-1)*gm*pr(i,j,k)/ro(i,j,k))
      enddo
      enddo
      enddo

      call cipdxsrc(rodxh,rodyh,rodzh,ro,roh,vx,vy,vz
     &                   ,dt,dx,dy,dz,ix,jx,kx)
      call cipdxsrc(tedxh,tedyh,tedzh,te,teh,vx,vy,vz
     &                   ,dt,dx,dy,dz,ix,jx,kx)
      call cipdxsrc(vxdxmh,vxdymh,vxdzmh,vxm,vxmh,vxm,vynx,vznx
     &                    ,dt,dx,dy,dz,ix,jx,kx)
      call cipdxsrc(vydxmh,vydymh,vydzmh,vym,vymh,vxny,vym,vzny
     &                    ,dt,dx,dy,dz,ix,jx,kx)
      call cipdxsrc(vzdxmh,vzdymh,vzdzmh,vzm,vzmh,vxnz,vynz,vzm
     &                    ,dt,dx,dy,dz,ix,jx,kx)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      call cipadv(roh,rodxh,rodyh,rodzh,vx,vy,vz,0,0,0
     &               ,dt,dxm,dym,dzm,ix,jx,kx)
      call cipadv(teh,tedxh,tedyh,tedzh,vx,vy,vz,0,0,0
     &               ,dt,dxm,dym,dzm,ix,jx,kx)
      call cipadv(vxmh,vxdxmh,vxdymh,vxdzmh,vxm,vynx,vznx,1,0,0
     &                ,dt,dx,dym,dzm,ix,jx,kx)
      call cipadv(vymh,vydxmh,vydymh,vydzmh,vxny,vym,vzny,0,1,0
     &                ,dt,dxm,dy,dzm,ix,jx,kx)
      call cipadv(vzmh,vzdxmh,vzdymh,vzdzmh,vxnz,vynz,vzm,0,0,1
     &                ,dt,dxm,dym,dz,ix,jx,kx)

c----------------------------------------------------------------------|
c--- MOC/CT Method
c----------------------------------------------------------------------|
c--- ### important ### Use (n+1)-step values for vx , vy & vz !

      call mocx(vxcy,vxcz,bxcy,bxcz,ro,vxmh,bxm,vymh,bym,vzmh,bzm
     &              ,dt,dxm,dym,dzm,ix,jx,kx)
      call mocy(vycx,vycz,bycx,bycz,ro,vxmh,bxm,vymh,bym,vzmh,bzm
     &              ,dt,dxm,dym,dzm,ix,jx,kx)
      call mocz(vzcx,vzcy,bzcx,bzcy,ro,vxmh,bxm,vymh,bym,vzmh,bzm
     &              ,dt,dxm,dym,dzm,ix,jx,kx)
      call ctransptx(bxmh,vxcy,vxcz,bxcy,bxcz,vycx,vzcx,bycx,bzcx
     &                   ,dt,dy,dz,ix,jx,kx)
      call ctranspty(bymh,vycx,vycz,bycx,bycz,vxcy,vzcy,bxcy,bzcy
     &                   ,dt,dx,dz,ix,jx,kx)
      call ctransptz(bzmh,vzcx,vzcy,bzcx,bzcy,vxcz,bxcz,vycz,bycz
     &                   ,dt,dx,dy,ix,jx,kx)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
        ro(i,j,k)=roh(i,j,k)
        rodx(i,j,k)=rodxh(i,j,k)
        rody(i,j,k)=rodyh(i,j,k)
        rodz(i,j,k)=rodzh(i,j,k)
        te(i,j,k)=teh(i,j,k)
        tedx(i,j,k)=tedxh(i,j,k)
        tedy(i,j,k)=tedyh(i,j,k)
        tedz(i,j,k)=tedzh(i,j,k)
        vxm(i,j,k)=vxmh(i,j,k)
        vxdxm(i,j,k)=vxdxmh(i,j,k)
        vxdym(i,j,k)=vxdymh(i,j,k)
        vxdzm(i,j,k)=vxdzmh(i,j,k)
        vym(i,j,k)=vymh(i,j,k)
        vydxm(i,j,k)=vydxmh(i,j,k)
        vydym(i,j,k)=vydymh(i,j,k)
        vydzm(i,j,k)=vydzmh(i,j,k)
        vzm(i,j,k)=vzmh(i,j,k)
        vzdxm(i,j,k)=vzdxmh(i,j,k)
        vzdym(i,j,k)=vzdymh(i,j,k)
        vzdzm(i,j,k)=vzdzmh(i,j,k)
        bxm(i,j,k)=bxmh(i,j,k)
        bym(i,j,k)=bymh(i,j,k)
        bzm(i,j,k)=bzmh(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
        pr(i,j,k) = ro(i,j,k)*te(i,j,k)/gm
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=2,ix
        vx(i,j,k) = (vxm(i-1,j,k)+vxm(i,j,k))/2
        bx(i,j,k) = (bxm(i-1,j,k)+bxm(i,j,k))/2
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx
      do i=1,ix
        vy(i,j,k) = (vym(i,j-1,k)+vym(i,j,k))/2
        by(i,j,k) = (bym(i,j-1,k)+bym(i,j,k))/2
      enddo
      enddo
      enddo
      do k=2,kx
      do j=1,jx
      do i=1,ix
        vz(i,j,k) = (vzm(i,j,k-1)+vzm(i,j,k))/2
        bz(i,j,k) = (bzm(i,j,k-1)+bzm(i,j,k))/2
      enddo
      enddo
      enddo

      return
      end
