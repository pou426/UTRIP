c======================================================================|
      subroutine cfl_h_d(dt,safety,dtmin,merr,gm,ro,pr,vx,vy
     &   ,x,dx,ix,dy,jx)
c======================================================================|
c 
c NAME  cfl_h_s
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * hydrodynamic equations
c        * Spherical coordinate, axis-symmetry
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
c    dx(ix),dy(jx): [double] grid spacing
c    x(ix): [double] coordinate
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension x(ix)
      dimension ro(ix,jx), pr(ix,jx), vx(ix,jx), vy(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|

c     dtmin=2.0e-10
c     safety = 0.40

      dt=1.e20
      imin   = 0
      jmin   = 0

      do j=1,jx
      do i=1,ix
         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)
         cs2 = gm*pr(i,j)/ro(i,j)
         dtcfl = min(dx(i),x(i)*dy(j))/sqrt(v2+cs2)
         dtq(i,j)=safety*dtcfl
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
