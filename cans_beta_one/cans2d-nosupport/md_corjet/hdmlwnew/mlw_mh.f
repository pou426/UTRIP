c======================================================================|
      subroutine mlw_mh(rof,eef,rxf,ryf,bxf,byf
     &  ,ro2,ee2,rx2,ry2,bx2,by2,az2
     &  ,roh,eeh,rxh,ryh,bxh,byh,azh
     &  ,gm,dt,dx,dxm,dy,dym,ix,jx,mstage)
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

      dimension rof(ix,jx),eef(ix,jx),rxf(ix,jx),ryf(ix,jx)
     &         ,bxf(ix,jx),byf(ix,jx)
      dimension prf(ix,jx),vxf(ix,jx),vyf(ix,jx),ezf(ix,jx)

      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
     &         ,bxh(ix,jx),byh(ix,jx),azh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx),ezh(ix,jx)

      dimension ro2(ix,jx),ee2(ix,jx),rx2(ix,jx),ry2(ix,jx)
     &         ,bx2(ix,jx),by2(ix,jx),az2(ix,jx)

      dimension fx(ix,jx)
      dimension fy(ix,jx)
      dimension ss(ix,jx)

      double precision mu

c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      mu=4.d0*pi

c----------------------------------------------------------------------|
c     calculate energy from pressure
c----------------------------------------------------------------------|
      call rv2vv(rof,rxf,ryf,vxf,vyf,ix,jx,+1)
      call ee2pr(prf,eef,rof,rxf,ryf,bxf,byf,gm,ix,jx,+1)

      call bb2ez(ezf,vxf,vyf,bxf,byf,ix,jx)
c----------------------------------------------------------------------|
c     step intermediate results for flux calculation
c----------------------------------------------------------------------|
c---  density ---
      call getfro(fx,fy,rof,vxf,vyf,ix,jx)
      if (mstage.eq.1) then
        call mlwh(roh,ro2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(ro2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

c---  energy ---
      call getfee(fx,fy,rof,prf,vxf,vyf,bxf,byf,ezf,gm,mu,ix,jx)
      if (mstage.eq.1) then
        call mlwh(eeh,ee2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(ee2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

c---  x-momentum ---
      call getfrx(fx,fy,rof,prf,vxf,vyf,bxf,byf,mu,ix,jx)
      if (mstage.eq.1) then
        call mlwh(rxh,rx2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(rx2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

c---  y-momentum ---
      call getfry(fx,fy,rof,prf,vxf,vyf,bxf,byf,mu,ix,jx)
      if (mstage.eq.1) then
        call mlwh(ryh,ry2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      else
        call mlwf(ry2,dt,fx,dx,dxm,fy,dy,dym,ix,jx)
      endif

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
