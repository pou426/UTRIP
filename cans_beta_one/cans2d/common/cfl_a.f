c======================================================================|
      subroutine cfl_a(dt,safety,dtmin,merr,vx,vy,dx,ix,dy,jx)
c======================================================================|
c 
c NAME  cfl_a
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * advection equations
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c
c INPUTS
c    ix,jx: [integer] dimension size
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    dx(ix),dy(jx): [double] grid spacing
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension vx(ix,jx), vy(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|

c     dtmin=2.0e-10
c     safety = 0.40
c     safety = 0.20

      dt=1.e20
      imin   = 0
      jmin   = 0

      do j=1,jx
      do i=1,ix
         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)
         dtcfl = min(dx(i),dy(j))/sqrt(v2)
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
