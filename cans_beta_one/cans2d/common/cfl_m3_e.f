c======================================================================|
      subroutine cfl_m3_e(dt,safety,dtmin,merr,gm
     &   ,ro,pr,vx,vy,vz,bx,by,bz,et,dx,ix,dy,jx)
c======================================================================|
c 
c NAME  cfl_m3_e
c
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * 3-component MHD equations
c        * resistivity
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
c    vz(ix,jx): [double] velocity
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    bz(ix,jx): [double] magnetic field
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
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx),vz(ix,jx)
     &         ,bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension et(ix,jx)
      dimension dtq(ix,jx)
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi4i=1./pi/4.

c     dtmin=2.0e-10
c     safety = 0.40

      etmin=1.e-5
      dt=1.e20
      imin   = 0
      do j=1,jx
      do i=1,ix
         dmesh=min(dx(i),dy(j))

         b2 = bx(i,j)*bx(i,j)+by(i,j)*by(i,j)+bz(i,j)*bz(i,j)
         ca2 = b2*pi4i/ro(i,j)
         cs2 = gm*pr(i,j)/ro(i,j)
         v2 = vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j)
         dtcfl = dmesh/sqrt(v2+cs2+ca2)
         dtwav=safety*dtcfl

         dtdif = safety*0.5*dmesh**2/max(et(i,j),etmin)
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
