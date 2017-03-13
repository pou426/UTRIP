c======================================================================|
      subroutine ee2eh(eh,ee,ro,vx,vy,vz,bx,by,bz,mu,ix,jx,kx)
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

      dimension eh(ix,jx,kx)
      dimension ro(ix,jx,kx),ee(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)

      double precision mu
c----------------------------------------------------------------------|
c     energy -> en-thermal

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ek= (vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)*ro(i,j,k)/2.d0
         em= (bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)/(2.d0*mu)
         eh(i,j,k)=ee(i,j,k)-ek-em
      enddo
      enddo
      enddo


      return
      end
