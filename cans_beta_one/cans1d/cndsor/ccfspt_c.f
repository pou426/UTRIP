c======================================================================|
      subroutine ccfspt_c(cmat,src,rkap0,gm,dt,te,ro,sc,scm,dx,dxm,ix)
c======================================================================|
c
c NAME  ccfspt_c
c
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Spitzer type
c        * non-uniform cross section
c
c OUTPUTS
c    cmat(ix,3): [double] coefficient matrix of heat conduction eq
c    src(ix): [double] source vector of heat conduction eq
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ix: [integer] dimension size
c    te(ix): [double] temperature
c    ro(ix): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    sc(ix), scm(ix) : [double] cross section
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension cmat(ix,3),src(ix)

      dimension ro(ix)
      dimension te(ix)

      dimension a0t(ix),axx(ix)

      dimension dx(ix),dxm(ix)

      dimension sc(ix),scm(ix)
c----------------------------------------------------------------------|

      telim=1.e4

      do i=1,ix-1
c        te00=0.5*(te(i)+te(i+1))
         te00=sqrt(te(i)*te(i+1))
         axx(i)=rkap0*scm(i)*sqrt(min(te00,telim))**5
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
