c======================================================================|
      subroutine model(ro,vx,cs2,g0,sc,scm,dv,dvm,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),vx(ix)
      dimension sc(ix),scm(ix),dv(ix),dvm(ix)
      dimension gr2(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      cs2=1.0d0

      g0=1.0d0

      xmin=1.d-3
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------

      dx0=1.d-6
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      i1=izero+10
      apx1=1.05d0
      dxmax1=0.1d0

      do i=izero+i1+1,ix
         dxm(i)=apx1*dxm(i-1)
         if(dxm(i).gt.dxmax1) dxm(i)=dxmax1
      enddo

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
        scm(i)=(x(i)+0.5d0*dxm(i))**2
      enddo

      do i=1,ix
         dv(i) = sc(i)
         dvm(i) = scm(i)
      enddo
c     dv(izero  )=1.d0/3.d0*((dxm(izero-1)+dxm(izero))/2)**2
c     dv(izero+1)=1.d0/3.d0*((dxm(izero+1)+dxm(izero))/2)**2
c     dvm(izero )=2.d0/3.d0*(dxm(izero)/2)**2

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      p4=pi*4.

      gr2(izero)=g0*x(izero)**2
      ro(izero)=1.d0

      do i=izero+1,ix
        rom=ro(i-1)
        gr2m=gr2(i-1)
        xm=x(i-1)
        rk01=dxm(i-1)*(-rom*gr2m/xm**2)
        rk02=dxm(i-1)*(p4*rom*xm**2)

        rom=ro(i-1)+dxm(i-1)/2*rk01
        gr2m=gr2(i-1)+dxm(i-1)/2*rk02
        xm=x(i-1)+dxm(i-1)/2
        rk11=dxm(i-1)*(-rom*gr2m/xm**2)
        rk12=dxm(i-1)*(p4*rom*xm**2)

        rom=ro(i-1)+dxm(i-1)/2*rk11
        gr2m=gr2(i-1)+dxm(i-1)/2*rk12
        xm=x(i)+dxm(i-1)/2
        rk21=dxm(i-1)*(-rom*gr2m/xm**2)
        rk22=dxm(i-1)*(p4*rom*xm**2)

        rom=ro(i-1)+dxm(i-1)*rk21
        gr2m=gr2(i-1)+dxm(i-1)*rk22
        xm=x(i-1)+dxm(i-1)
        rk31=dxm(i-1)*(-rom*gr2m/xm**2)
        rk32=dxm(i-1)*(p4*rom*xm**2)

        ro(i) =ro(i-1)+1./6.*(rk01+2.*rk11+2.*rk21+rk31)
        gr2(i)=gr2(i-1)+1./6.*(rk02+2.*rk12+2.*rk22+rk32)
      enddo
      do i=izero-1,1,-1
        ro(i)=ro(izero)
      enddo

      amp=4.d0
      do i=1,ix
         ro(i) = amp*ro(i)
         vx(i) = 0.0d0
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'amp',amp)

      
      return
      end
