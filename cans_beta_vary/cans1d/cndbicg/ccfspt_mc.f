c======================================================================|
      subroutine ccfspt_mc(cmat,src,rkap0,gm,dt,te,ro,bxm,by,sc,scm
     &  ,dx,dxm,ix)
c======================================================================|
c
c NAME  ccfspt_mc
c
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Spitzer type
c        * magnetic field 2 components
c        * non-uniform cross section
c
c OUTPUTS
c    cmat(ix,3): [double] coefficient matrix of heat conduction eq
c    src(ix): [double] source vector of heat conduction eq
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    te(ix): [double] temperature
c    ro(ix): [double] density
c    bxm(ix) : [double] magnetic field
c    by(ix) : [double] magnetic field
c    sc(ix), scm(ix) : [double] cross section
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,3),src(ix)
      dimension ro(ix)
      dimension te(ix)
      dimension by(ix),bxm(ix)
      dimension sc(ix),scm(ix)
      dimension a0t(ix),axx(ix)
      dimension dx(ix),dxm(ix)
c----------------------------------------------------------------------|
      telim=1.e4

      do i=1,ix-1
c        te00=0.5*(te(i)+te(i+1))
         te00=sqrt(te(i)*te(i+1))
         bt00=0.5*(by(i)+by(i+1))
         ebb=bxm(i)**2/(bt00**2+bxm(i)**2)
         axx(i)=rkap0*scm(i)*sqrt(min(te00,telim))**5*ebb
      enddo
      do i=1,ix
         a0t(i)=sc(i)*ro(i)/(gm-1.)
      enddo

       do i=2,ix-1
          aa0=a0t(i)/dt
          src(i)=aa0*te(i)
       enddo
       do i=2,ix-1
          aa0=a0t(i)/dt
          aaxp=axx(i)/dxm(i)/dx(i)
          aaxm=axx(i-1)/dxm(i-1)/dx(i)
          cmat(i,3)=-aaxp
          cmat(i,2)=-aaxm
          cmat(i,1)=aa0+aaxp+aaxm
       enddo


      return
      end
