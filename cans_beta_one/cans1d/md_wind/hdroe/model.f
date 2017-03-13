c======================================================================|
      subroutine model(ro,pr,vx,gm,gx,gxm,sc,scm,dv,margin,x,ix
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix)
      dimension dxm(ix)

      dimension ro(ix),pr(ix),vx(ix)

      dimension sc(ix),scm(ix),dv(ix)
      dimension gx(ix),gxm(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.05d0
      rstar=10.d0
      rmax=300.d0
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------
      dx0=(rmax-rstar)/real(ix-margin*2)
      ratx=1.001d0

c-----------------------------------------------------------------------
c      dxm,x

      do i=1,ix
         dxm(i)=dx0
      enddo
      do i=2,ix
         dxm(i)=dxm(i-1)*ratx
      enddo

      izero=margin
      x(izero)=rstar-dxm(izero)/2

        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo

c----------------------------------------------------------------------|
c       cross section of flux tube
c----------------------------------------------------------------------|

      do i=1,ix
        sc(i)=x(i)**2
        scm(i)=(x(i)+dxm(i)/2)**2
      enddo

      do i=1,ix
         dv(i) = sc(i)
      enddo

c----------------------------------------------------------------------|
c       gravity
c----------------------------------------------------------------------|
      do i=1,ix
        gx(i)=-1/gm*rstar**2/x(i)**2
        gxm(i)=-1/gm*rstar**2/(x(i)+dxm(i)/2)**2
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      rhs=2*rstar
      roism=1.d-8
      teism=1.d-2

      do i=1,ix
        vx(i)  = 0.d0
        ro(i)  = exp(rstar*(rstar/x(i)-1))*(1-tanh((x(i)-rhs)/0.5d0))/2
     &            + roism*(1+tanh((x(i)-rhs)/0.5d0))/2
c       te = 1+(teism-1)*(1+tanh((x(i)-rhs)/0.5d0))/2
        te = (teism+(1-teism)*(1-(x(i)-rhs)/(rmax-rhs)))
c    &         *(1+tanh((x(i)-rhs)/0.5d0))/2
        if (x(i).le.rhs) te = 1.d0
        pr(i)  = ro(i)*te/gm
        if (ro(i).le.1.d-14) ro(i)=1.d-14 
        if (pr(i).le.1.d-14) pr(i)=1.d-14 
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'rstar',rstar)
      call dacputparamd(mf_params,'rhs',rhs)
      call dacputparamd(mf_params,'roism',roism)
      call dacputparamd(mf_params,'teism',teism)


      return
      end
