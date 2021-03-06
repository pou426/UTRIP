c======================================================================|
      subroutine cfl_mt_ce(dt,safety,dtmin,merr,cs2,ro,vx,vy,vz,bx,by,bz
     &         ,et,x,dx,ix,dy,jx,dz,kx)
c======================================================================|
c 
c NAME  cfl_m_ce
c 
c PURPOSE
c    determine time step such that it satisfies CFL condition.
c        * isothermal MHD equations 
c        * Cylindrical coordinate
c        * Resistivity 
c        
c OUTPUTS
c    dt: [double] delta time
c    merr: [integer] error code, merr=0 is nominal.
c    
c INPUTS
c    ix,jx,kx: [integer] dimension size
c    ro(ix,jx,kx): [double] density
c    vx(ix,jx,kx): [double] velocity
c    vy(ix,jx,kx): [double] velocity 
c    bx(ix,jx,kx): [double] magnetic field
c    by(ix,jx,kx): [double] magnetic field
c    et(ix,jx,kx): [double] resistivity
c    dx(ix),dy(jx): [double] grid spacing
c    x(ix): [double] coordinate
c    cs2: [double] square of sound speed
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx),dz(kx)
      dimension x(ix)
      dimension ro(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension et(ix,jx,kx)
      dimension dtq(ix,jx,kx)
c----------------------------------------------------------------------|

      pai = acos(-1.0)
      onep4=0.25/pai

c     dtmin=2.0e-10
c     safety = 0.40

      dt=1.e20
      imin   = 0
      jmin   = 0
      kmin   = 0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         dmesh = min(dx(i),x(i)*dy(j),dz(k))

         b2 = bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         ca2 = b2*onep4/ro(i,j,k)
         v2 = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         dtcfl = dmesh/sqrt(v2+cs2+ca2)
         dtwav=safety*dtcfl

         dtdif = safety*0.5*dmesh**2/max(et(i,j,k),1.d-10)
         dtq(i,j,k)=min(dtwav,dtdif)

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
