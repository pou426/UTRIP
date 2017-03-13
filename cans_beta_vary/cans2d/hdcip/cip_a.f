c======================================================================|
      subroutine cip_a(ro,rodx,rody,dt,vx,vxm,vy,vym
     &             ,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  cip_a
c
c PURPOSE
c    solve eqs. by modified Roe + MUSCL-TVD  method with effects of
c        * simple advection
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
      dimension vx(ix,jx),vxm(ix,jx),vy(ix,jx),vym(ix,jx)
      dimension ro(ix,jx),rodx(ix,jx),rody(ix,jx)
      dimension roh(ix,jx),rodxh(ix,jx),rodyh(ix,jx)

c----------------------------------------------------------------------|
c--- preparation
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
        roh(i,j)=ro(i,j)
        rodxh(i,j)=rodx(i,j)
        rodyh(i,j)=rody(i,j)
      enddo
      enddo

c----------------------------------------------------------------------|
c--- non advection & source term ---
c----------------------------------------------------------------------|
      do j=2,jx
      do i=2,ix
        dvx=vxm(i,j)-vxm(i-1,j)
        dvy=vym(i,j)-vym(i-1,j)
        divv=dvx/dx(i)+dvy/dy(j)
        roh(i,j)=roh(i,j)+dt*(-divv*ro(i,j))
      enddo
      enddo

      call cipdxsrc(rodxh,rodyh,ro,roh,vx,vy,dt,dx,dy,ix,jx)

c----------------------------------------------------------------------|
c--- advection phase ---
c----------------------------------------------------------------------|

      call cipadv(roh,rodxh,rodyh,vx,vy,0,0,dt,dxm,dym,ix,jx)

c----------------------------------------------------------------------|
c--- ending
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
        ro(i,j)=roh(i,j)
        rodx(i,j)=rodxh(i,j)
        rody(i,j)=rodyh(i,j)
      enddo
      enddo


      return
      end
