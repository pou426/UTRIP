c======================================================================|
      subroutine cip_h(ro,pr,vx,vxm,vy,vym,te
     &             ,rodx,rody,tedx,tedy,vxdxm,vxdym,vydxm,vydym
     &             ,dt,cvis,gm,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  cip_h
c
c PURPOSE
c    solve eqs. by CIP method
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
      dimension rodx(ix,jx),rody(ix,jx),tedx(ix,jx),tedy(ix,jx)
      dimension vxdxm(ix,jx),vxdym(ix,jx),vydxm(ix,jx),vydym(ix,jx)
      dimension roh(ix,jx),teh(ix,jx),vxmh(ix,jx),vymh(ix,jx)
      dimension rodxh(ix,jx),rodyh(ix,jx),tedxh(ix,jx),tedyh(ix,jx)
      dimension vxdxmh(ix,jx),vxdymh(ix,jx),vydxmh(ix,jx),vydymh(ix,jx)
      dimension pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension vxn(ix,jx),vyn(ix,jx)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

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
      enddo
      enddo
      do j=2,jx
      do i=1,ix
        vy(i,j) = (vym(i,j-1)+vym(i,j))/2
      enddo
      enddo
      do j=1,jx-1
      do i=2,ix
        vxn(i,j) = (vxm(i-1,j)+vxm(i,j)+vxm(i-1,j+1)+vxm(i,j+1))/4
      enddo
      enddo
      do j=2,jx
      do i=1,ix-1
        vyn(i,j) = (vym(i,j-1)+vym(i,j)+vym(i+1,j-1)+vym(i+1,j))/4
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
      do j=1,jx
      do i=1,ix-1
        vxmh(i,j)=vxmh(i,j)
     &       +dt*(-(pr(i+1,j)-pr(i,j))/dxm(i)/((ro(i,j)+ro(i+1,j))/2))
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
        vymh(i,j)=vymh(i,j)
     &       +dt*(-(pr(i,j+1)-pr(i,j))/dym(j)/((ro(i,j)+ro(i,j+1))/2))
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
      call cipdxsrc(vxdxmh,vxdymh,vxm,vxmh,vxm,vym,dt,dx,dy,ix,jx)
      call cipdxsrc(vydxmh,vydymh,vym,vymh,vxm,vym,dt,dx,dy,ix,jx)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      call cipadv(roh,rodxh,rodyh,vx,vy,0,0,dt,dxm,dym,ix,jx)
      call cipadv(teh,tedxh,tedyh,vx,vy,0,0,dt,dxm,dym,ix,jx)
      call cipadv(vxmh,vxdxmh,vxdymh,vxm,vyn,1,0,dt,dx,dym,ix,jx)
      call cipadv(vymh,vydxmh,vydymh,vxn,vym,0,1,dt,dxm,dy,ix,jx)

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
      enddo
      enddo
      do j=2,jx
      do i=1,ix
        vy(i,j) = (vym(i,j-1)+vym(i,j))/2
      enddo
      enddo


      return
      end
