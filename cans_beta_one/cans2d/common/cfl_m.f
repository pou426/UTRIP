c======================================================================|
      subroutine cfl_m(dt,safety,dtmin,merr,gm,ro,pr,vx,vy,bx,by
     &   ,dx,ix,dy,jx)
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
c    ix,jx: [integer] dimension size
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    dx(ix),dy(jx): [double] grid spacing
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix)
      dimension dy(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
     &         ,bx(ix,jx),by(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c     dtmin=2.0e-10
c     safety = 0.40

      dt=1.e20
      imin   = 0
      do j=1,jx
      do i=1,ix
         b2 = bx(i,j)*bx(i,j)+by(i,j)*by(i,j)
         ca2 = b2*pi4i/ro(i,j)
         cs2 = gm*pr(i,j)/ro(i,j)
         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)
         dtcfl = min(dx(i),dy(j))/sqrt(v2+cs2+ca2)
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
