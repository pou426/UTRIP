c======================================================================|
      subroutine cfl_m(dt,safety,dtmin,merr,gm,bx,ro,pr,vx,vy,by,dx,ix)
c======================================================================|
c 
c NAME  cfl_m
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * MHD equations
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c
c INPUTS
c    ix: [integer] dimension size
c    ro(ix): [double] density
c    pr(ix): [double] pressure
c    vx(ix): [double] velocity
c    vy(ix): [double] velocity
c    bx(ix): [double] magnetic field
c    by(ix): [double] magnetic field
c    dx(ix): [double] grid spacing
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
c     determine time step such that it satisfies cfl condition.
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension dx(ix)
      dimension ro(ix)
      dimension pr(ix)
      dimension vx(ix)
      dimension vy(ix)
      dimension bx(ix)
      dimension by(ix)
      dimension dtq(ix)
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=0.25/pi
c----------------------------------------------------------------------|
      dt=1.e20
      imin   = 0
      do i=2,ix-1
         onero = 1.0/ro(i)
         ca2 = (bx(i)*bx(i)+by(i)*by(i))*pi4i*onero
         cs2 = gm*pr(i)*onero
         dtcfl = dx(i)/(abs(vx(i))+sqrt(cs2+ca2))
         dtq(i)=safety*dtcfl
      enddo

      do i=2,ix-1
        if(dtq(i).lt.dt) then
          imin=i
          dt=dtq(i)
        endif
      enddo
c----------------------------------------------------------------------|
c     write the point where dt is smaller than critical value    
c----------------------------------------------------------------------|
      merr=0

      if (dt.lt.dtmin) then
         merr=9001
         write(6,*) '  ### stop due to small dt, less than dtmin ###'
         write(6,620) dt,dtmin,imin
 620     format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5) 
      endif

      return
      end
