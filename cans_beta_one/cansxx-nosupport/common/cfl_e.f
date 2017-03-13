c======================================================================|
      subroutine cfl_e(dt,safety,dtmin,et,dx,dy,ix,jx,merr)
c======================================================================|
c 
c NAME  cfl_m_e
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * MHD equations
c        * Resistivity
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
c    et(ix,jx): [double] resistivity
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
      dimension et(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|

      etmin=1.d-20
      dt=1.d20

      do j=1,jx
      do i=1,ix
         dmesh=min(dx(i),dy(j))
         dtq(i,j) = safety*0.5d0*dmesh**2/max(et(i,j),etmin)
      enddo
      enddo

      imin   = 0
      jmin   = 0
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
