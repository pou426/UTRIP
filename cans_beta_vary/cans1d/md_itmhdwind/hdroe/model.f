c======================================================================|
      subroutine model(ro,vx,vy,by
     &   ,vstar,bx,bxm,cs2,gx,gxm,rr,rrm,sc,scm,dv,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),vx(ix),vy(ix),by(ix)

      dimension sc(ix),scm(ix),dv(ix)
      dimension bx(ix),bxm(ix)
      dimension gx(ix),gxm(ix)
      dimension rr(ix),rrm(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      cs2=1.d0
      rstar=6.d0
      vstar=0.15d0
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------

      dx0=10.d0/real(ix-margin*2)
      ratx=1.006d0

c-----------------------------------------------------------------------
c      dxm,x

      do i=1,ix
         dxm(i)=dx0
      enddo
      do i=2,ix
         dxm(i)=dxm(i-1)*ratx
      enddo

      izero=margin+1
      x(izero)=rstar+dxm(izero)/2

        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo

c----------------------------------------------------------------------|
c       gravity
c----------------------------------------------------------------------|
      do i=1,ix
        gx(i)=-rstar**2/x(i)**2
        gxm(i)=-rstar**2/(x(i)+0.5*dxm(i))**2
      enddo
c----------------------------------------------------------------------|
c   curvature of magnetic field
c----------------------------------------------------------------------|
      do i=1,ix
        rr(i)=x(i)
        rrm(i)=x(i)+dxm(i)/2
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      rhs=2.0d0*rstar
      roism=1.d-8
      teism=1.d-2

      do i=1,ix
         vx(i)  = 0.d0
         ro(i)  = exp(rstar*(rstar/x(i)-1))*(1-tanh((x(i)-rhs)/0.5))/2
     &            + roism*(1+tanh((x(i)-rhs)/0.5))/2
      enddo
c----------------------------------------------------------------------|
c       magnetic field
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      betai=0.1d0
      bx0=sqrt(betai*cs2*ro(izero)*8*pi)
      do i=1,ix
        vy(i)=0.d0
        bx(i)=bx0*x(1)**2/x(i)**2
        bxm(i)=bx0*x(1)**2/(x(i)+dxm(i)/2)**2
        by(i)=0.d0
      enddo
c----------------------------------------------------------------------|
c       cross section of flux tube
c----------------------------------------------------------------------|

      do i=1,ix
        sc(i)=1.0d0/bx(i)
        scm(i)=1.0d0/bxm(i)
      enddo
      do i=1,ix
        dv(i)=sc(i)
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'rstar',rstar)
      call dacputparamd(mf_params,'rhs',rhs)
      call dacputparamd(mf_params,'roism',roism)
      call dacputparamd(mf_params,'teism',teism)
      call dacputparamd(mf_params,'betai',betai)


      return
      end
