c======================================================================|
      subroutine cfl_h_kv(dt,merr,gm,ro,pr,vx,vy,rkap,visc,dx,ix,dy,jx)
c======================================================================|
c 
c NAME  cfl_h_kv
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * HD equations
c        * Thermal conduction
c        * Viscosity
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c
c INPUTS
c    ix,jx: [integer] dimension size
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    rkap     : [double] thermal conductivity
c    visc     : [double] viscosity
c    dx(ix),dy(jx): [double] grid spacing
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2003-5-1 H. Isobe
c
c----------------------------------------------------------------------|
      implicit real*8 (a-h,o-z)
      dimension dx(ix)
      dimension dy(jx)
      dimension ro(ix,jx), pr(ix,jx), vx(ix,jx), vy(ix,jx)
      dimension rkap(ix,jx),visc(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|

      dtmin=2.0e-10
      safety = 0.10

      dt=1.e20
      imin   = 0
      jmin   = 0

      do j=1,jx
      do i=1,ix
         dmesh=min(dx(i),dy(j))

         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)
         cs2 = gm*pr(i,j)/ro(i,j)
         dtcfl = dmesh/sqrt(v2+cs2)
         dtwav=safety*dtcfl
         
         dtdif = safety*0.5*dmesh**2/max(rkap(i,j),visc(i,j))
         dtq(i,j)=min(dtwav,dtdif)
      enddo
      enddo

      do j=1,jx
      do i=1,ix
        if(dtq(i,j).lt.dt) then
          imin=i
          jmin=j
          dt=dtq(i,j)
        endif
      enddo
      enddo
c----------------------------------------------------------------------|
c     write the point where dt is smaller than critical value    
c----------------------------------------------------------------------|
      merr=0

      if (dt.lt.dtmin) then
         merr=9001
         write(6,*) '  ### stop due to small dt, less than dtmin ###'
         write(6,620) dt,dtmin,imin,jmin
 620     format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5
     &          ,' j = ',i5) 
      endif

      return
      end





