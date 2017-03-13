c======================================================================|
      subroutine mlw_g(rof,eef,rxf,gxf,ee2,rx2,eeh,rxh
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

      dimension rof(ix,jx),eef(ix,jx),rxf(ix,jx)
      dimension eeh(ix,jx),rxh(ix,jx)
      dimension ee2(ix,jx),rx2(ix,jx)

      dimension ss(ix,jx)

      dimension gxf(ix,jx)

c----------------------------------------------------------------------|
c     gravity
c----------------------------------------------------------------------|
c---  x-momentum ---
      do j=1,jx
      do i=1,ix       
         ss(i,j)= rof(i,j)*gxf(i,j)
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwsh(rxh,rx2,dt,ss,ix,jx)
      else
        call mlwsf(rx2,dt,ss,dx,dxm,dy,dym,ix,jx)
      endif

c---  x-comp. energy ---
      do j=1,jx
      do i=1,ix       
         ss(i,j)= rxf(i,j)*gxf(i,j)
      enddo
      enddo
      if (mstage.eq.1) then
        call mlwsh(eeh,ee2,dt,ss,ix,jx)
      else
        call mlwsf(ee2,dt,ss,dx,dxm,dy,dym,ix,jx)
      endif


      return
      end
