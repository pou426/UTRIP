c======================================================================|
      subroutine model(ro,pr,vx,gm,sc,scm,dv,margin,x,ix
     &    ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix)
      dimension dxm(ix)

      dimension sc(ix),scm(ix),dv(ix)

      dimension ro(ix),pr(ix),vx(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      xmin=0.04d0
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------
      dx0=1.d0/dble(ix-margin*2)

      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      x(izero)=xmin+dxm(izero)/2.d0

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
        scm(i)=(x(i)+0.5*dxm(i))**2
      enddo

      do i=1,ix
         dv(i) = sc(i)
      enddo
      dv(izero  )=1.d0/3.d0*((dxm(izero-1)+dxm(izero))/2)**2
      dv(izero+1)=1.d0/3.d0*((dxm(izero+1)+dxm(izero))/2)**2

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      prism=1.d-8
      wexp=0.1d0

      do i=1,ix
         ro(i)  = 1.d0
         vx(i) = 0.0d0
         pr(i)  = prism+(1.-prism)*exp(-(x(i)/wexp)**2)
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'wexp',wexp)
      call dacputparamd(mf_params,'prism',prism)

      return
      end
