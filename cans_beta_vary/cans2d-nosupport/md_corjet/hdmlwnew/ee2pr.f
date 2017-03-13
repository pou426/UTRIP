c======================================================================|
      subroutine ee2pr(pr,ee,ro,rx,ry,bx,by,gm,ix,jx,mdir)
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

      dimension pr(ix,jx)
      dimension ro(ix,jx),ee(ix,jx),rx(ix,jx),ry(ix,jx)
     &         ,bx(ix,jx),by(ix,jx)

      double precision mu
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      mu = 4.d0*pi
c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
c     energy -> pressure
      if (mdir.gt.0) then

      do j=1,jx
      do i=1,ix
         rr= rx(i,j)**2+ry(i,j)**2
         bb= bx(i,j)**2+by(i,j)**2
         pr(i,j) = (gm-1.d0)*(ee(i,j) - 0.5d0*rr/ro(i,j) -bb/(2.d0*mu))
      enddo
      enddo

c     pressure -> energy
      else

      do j=1,jx
      do i=1,ix
         rr= rx(i,j)**2+ry(i,j)**2
         bb=bx(i,j)**2+by(i,j)**2
         ee(i,j) = pr(i,j)/(gm-1.d0)+0.5d0*rr/ro(i,j)+bb/(2.d0*mu)
      enddo
      enddo

      endif

      return
      end
