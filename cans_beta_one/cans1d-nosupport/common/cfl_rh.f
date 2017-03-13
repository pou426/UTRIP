c======================================================================|
      subroutine cfl_rh(dt,safety,dtmin,merr,gm,ro,pr,vx,hh,dx,ix)
c======================================================================|
c 
c NAME  cfl_rh
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * hydrodynamic equations
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
c    dx(ix): [double] grid spacing
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension dtq(ix)
      dimension hh(ix,4)
c----------------------------------------------------------------------|
      dt=1.e20
      imin   = 0
      do i=2,ix-1
         cs2 = gm*pr(i)/(ro(i)+gm/(gm-1)*pr(i))
         dtcfl = hh(i,2)*dx(i)/(abs(vx(i))+sqrt(cs2))
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
