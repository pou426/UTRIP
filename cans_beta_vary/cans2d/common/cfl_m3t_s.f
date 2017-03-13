c======================================================================|
      subroutine cfl_m3t_s(dt,safety,dtmin,merr,cs2
     &       ,ro,vx,vy,vz,bx,by,bz,x,dx,ix,dy,jx)
c======================================================================|
c 
c NAME  cfl_m3_s
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * isothermal 3-component MHD equations
c        * Spherical coordinate, axis-symmetry
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c
c INPUTS
c    ix,jx: [integer] dimension size
c    ro(ix,jx): [double] density
c    vx(ix,jx): [double] velocity
c    vy(ix,jx): [double] velocity
c    vz(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    bz(ix,jx): [double] magnetic field
c    dx(ix),dy(jx): [double] grid spacing
c    x(ix): [double] coordinate
c    cs2: [double] square of sound speed
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension x(ix) 
      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx),vz(ix,jx)
     &         ,bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c     dtmin=2.0e-10
c     dtmin=1.0e-5
c     safety = 0.40

      dt=1.e20
      imin   = 0
      do j=1,jx
      do i=1,ix
         b2 = bx(i,j)**2+by(i,j)**2+bz(i,j)**2
         ca2 = b2*pi4i/ro(i,j)
         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)**2
         dtcfl = min(dx(i),x(i)*dy(j))/sqrt(v2+cs2+ca2)
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
