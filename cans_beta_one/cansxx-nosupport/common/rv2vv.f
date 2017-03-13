c======================================================================|
      subroutine rv2vv(ro,rx,ry,rz,vx,vy,vz,ix,jx,kx,mdir)
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

      dimension ro(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)

c----------------------------------------------------------------------|
c     convert from total energy to pressure
c----------------------------------------------------------------------|
c     momentum -> velocity
      if (mdir.gt.0) then

      do k=1,kx
      do j=1,jx
      do i=1,ix
         vx(i,j,k)=rx(i,j,k)/ro(i,j,k)
         vy(i,j,k)=ry(i,j,k)/ro(i,j,k)
         vz(i,j,k)=rz(i,j,k)/ro(i,j,k)
      enddo
      enddo
      enddo

c     velocity -> momentum
      else

      do k=1,kx
      do j=1,jx
      do i=1,ix
         rx(i,j,k)=vx(i,j,k)*ro(i,j,k)
         ry(i,j,k)=vy(i,j,k)*ro(i,j,k)
         rz(i,j,k)=vz(i,j,k)*ro(i,j,k)
      enddo
      enddo
      enddo

      endif

      return
      end
