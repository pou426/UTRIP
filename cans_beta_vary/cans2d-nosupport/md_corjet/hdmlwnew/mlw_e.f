c======================================================================|
      subroutine mlw_e(bxf,byf,etf
     &  ,bx2,by2,az2,bxh,byh,azh
     &  ,dt,dx,dxm,dy,dym,ix,jx,mstage)
c======================================================================|
c
c NAME  mlw_m_g
c
c PURPOSE
c    solve eqs. by modified Lax-Wendroff method with effects of
c        * MHD
c        * gravity
c        * resistivity
c        
c INPUTS & OUTPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field 
c    az(ix,jx): [double] magnetic vector potential
c    
c OUTPUTS
c    None
c 
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c 
c    gx(ix,jx), gxm(ix,jx) : [double] gravity
c    gy(ix,jx), gym(ix,jx) : [double] gravity
c    et(ix), etm(ix): [double] resistivity
c    dx(ix), dxm(ix): [double] grid spacing
c    dy(jx), dym(jx): [double] grid spacing
c    dt: [double] delta time
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)

      dimension bxf(ix,jx),byf(ix,jx)

      dimension bxh(ix,jx),byh(ix,jx),azh(ix,jx)

      dimension bx2(ix,jx),by2(ix,jx),az2(ix,jx)

      dimension fx(ix,jx),fy(ix,jx),ss(ix,jx)

      dimension czf(ix,jx),ezf(ix,jx)
      dimension etf(ix,jx)
      double precision mu

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      mu=4.d0*pi

c----------------------------------------------------------------------|
c     resistivity
c----------------------------------------------------------------------|
      call bbtocz(czf,bxf,byf,dx,dy,ix,jx)

      do j=1,jx
      do i=1,ix
         ezf(i,j) = etf(i,j)*czf(i,j)
      enddo
      enddo

c---  x-magnetic ---
      call getfbx(fx,fy,ezf,ix,jx)
      if (mstage.eq.1) then
        call mlwh(bxh,bx2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(bx2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

c---  y-magnetic ---
      call getfby(fx,fy,ezf,ix,jx)
      if (mstage.eq.1) then
        call mlwh(byh,by2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(by2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

c---  z-magnetic potential ---
      call getsaz(ss,ezf,ix,jx)
      if (mstage.eq.1) then
        call mlwsh(azh,az2,dt,ss,ix,jx)
      else
        call mlwsf(az2,dt,ss,dx,dxm,dy,dym,ix,jx)
      endif


      return
      end
