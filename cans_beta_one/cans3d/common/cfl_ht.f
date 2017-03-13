c======================================================================|
      subroutine cfl_ht(dt,safety,dtmin,merr,cs2,vx,vy,vz
     &   ,dx,ix,dy,jx,dz,kx)
c======================================================================|
c
c NAME  cfl_ht
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * isothermal hydrodynamic equations
c
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c 
c INPUTS
c    ix,jx,kx: [integer] dimension size
c    vx(ix,jx,kx): [double] velocity
c    vy(ix,jx,kx): [double] velocity
c    vz(ix,jx,kx): [double] velocity
c    dx(ix),dy(jx): [double] grid spacing
c    cs2: [double] square of sound speed
c     
c HISTORY
c    written 2002-5-15 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx),dz(kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension dtq(ix,jx,kx)
c----------------------------------------------------------------------|

c     dtmin=2.0e-10
c     safety = 0.40

      dt=1.e20
      imin   = 0
      jmin   = 0
      kmin   = 0
      do k=1,kx
      do j=1,jx
      do i=1,ix
         v2 = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         dtcfl = min(dx(i),dy(j),dz(k))/sqrt(v2+cs2)
         dtq(i,j,k)=safety*dtcfl
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix
        if(dtq(i,j,k).lt.dt) then
          imin=i
          jmin=j
          kmin=k
          dt=dtq(i,j,k)
        endif
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     write the point where dt is smaller than critical value    
c----------------------------------------------------------------------|
      merr=0

      if (dt.lt.dtmin) then
         merr=9001
         write(6,*) '  ### stop due to small dt, less than dtmin ###'
         write(6,620) dt,dtmin,imin,jmin,kmin
 620     format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5
     &          ,' j = ',i5,' k = ',i5) 
      endif

      return
      end
