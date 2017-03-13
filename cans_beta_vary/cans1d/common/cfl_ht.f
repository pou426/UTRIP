c======================================================================|
      subroutine cfl_ht(dt,safety,dtmin,merr,cs2,vx,dx,ix)
c======================================================================|
c 
c NAME  cfl_ht
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * isotropic hydrodynamic equations
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c
c INPUTS
c    ix: [integer] dimension size
c    vx(ix): [double] velocity
c    dx(ix): [double] grid spacing
c    cs2: [double] square of sound speed
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix)
      dimension vx(ix)
      dimension dtq(ix)
c----------------------------------------------------------------------|
      dt=1.e20
      imin   = 0
      cs=sqrt(cs2)
      do i=2,ix-1
         v2 = vx(i)*vx(i)
         dtcfl = dx(i)/(abs(vx(i))+cs)
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
