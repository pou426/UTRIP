c======================================================================|
      subroutine model(ro,pr,vx,gm,rkap0,cool0,tecl0,decl0 
     &   ,gx,gxm,sc,scm,dv,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)

      dimension sc(ix),scm(ix),dv(ix)

      dimension ro(ix),pr(ix),vx(ix)

      dimension gx(ix),gxm(ix)
      dimension te1d(ix),pr1d(ix),ro1d(ix)

c----------------------------------------------------------------------|
c     parameters
c----------------------------------------------------------------------|
      ztr=12.5d0

      gm=5.d0/3.d0

      tenml=1.d4
      ronml=1.d17
      rlnml=2.d7

      call cndprm(rkap0,gm,tenml,ronml,rlnml)
      call htclprm(cool0,tecl0,decl0,gm,tenml,ronml,rlnml)
c-----------------------------------------------------------------------
c   grid
c-----------------------------------------------------------------------

      dx0=0.02d0
      i1=aint(1.3d0*ztr/dx0)
      apx1=1.02d0
      dxmax1=50.d0*dx0

      i0=0
      apx0=1.d0
      dxmax0=dx0


      do i=1,ix
         dxm(i)=dx0
      enddo

      if (i1.lt.ix) then
      do i=i1+1,ix
         dxm(i)=apx1*dxm(i-1)
         if(dxm(i).gt.dxmax1) dxm(i)=dxmax1
      enddo
      endif
      
      if (i0.gt.1) then
      do i=i0-1,1,-1
         dxm(i)=apx0*dxm(i+1)
         if(dxm(i).gt.dxmax0) dxm(i)=dxmax0
      enddo
      endif
      call bdsmpx(0,margin-1,dxm,ix)
       
      izero=margin
        x(izero)=-dxm(izero)/2
        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo

      xtop=(x(ix-margin)+dxm(ix-margin))

c----------------------------------------------------------------------|
c       cross section of flux tube
c----------------------------------------------------------------------|

      do i=1,ix
        sc(i)=1.d0
        scm(i)=1.d0
      enddo

      do i=1,ix
         dv(i) = sc(i)
      enddo

c----------------------------------------------------------------------|
c     parameters
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

      g0=1.0d0/gm
      ro0=1.0d0
      pr0=1.0d0/gm
      te0=1.0d0

      tpho=te0
      wtr=0.5d0
      tcor=200.0d0
      
c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      do i=1,ix
         gx(i)  = -g0*cos(0.5d0*pi*x(i)/xtop)
      enddo

      do i=1,ix
         gxm(i) = -g0*cos(0.5d0*pi*(x(i)+0.5d0*dxm(i))/xtop)
      enddo

      call bdspnx(0,margin,gx,ix)
      call bdsmnx(0,margin-1,gxm,ix)
      
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|
      do i=1,ix
            te1d(i)=tpho+0.5d0*(tcor-tpho)*(tanh((x(i)-ztr)/wtr)+1.)
      enddo

      i=izero+1
         ro1d(i) = (te1d(i-1)+0.5d0*gm*gxm(i-1)*dxm(i-1)/2)/
     &             (te1d(i)  -0.5d0*gm*gxm(i-1)*dxm(i-1)/2)*ro0

      do i=izero+2,ix
         ro1d(i) = (te1d(i-1)+0.5d0*gm*gxm(i-1)*dxm(i-1))/
     &             (te1d(i)  -0.5d0*gm*gxm(i-1)*dxm(i-1))*ro1d(i-1)
      enddo

      i=izero
         ro1d(i) = (te1d(i+1)-0.5d0*gm*gxm(i)*dxm(i)/2)/
     &             (te1d(i)  +0.5d0*gm*gxm(i)*dxm(i)/2)*ro0

      do i=izero-1,1,-1
         ro1d(i) = ro1d(i+1)
     &             *(te1d(i+1)-0.5d0*gm*gxm(i)*dxm(i))
     &             /(te1d(i)  +0.5d0*gm*gxm(i)*dxm(i))
      enddo

      do i=1,ix
         pr1d(i) = pr0*(ro1d(i)/ro0)*(te1d(i)/te0)
      enddo

c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      do i=1,ix
         ro(i)  = ro1d(i)
         vx(i)  = 0.0d0
         pr(i)  = pr1d(i)
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ztr',ztr)
      call dacputparamd(mf_params,'tcor',tcor)
      call dacputparamd(mf_params,'tenml',tenml)
      call dacputparamd(mf_params,'ronml',ronml)
      call dacputparamd(mf_params,'rlnml',rlnml)

      
      return
      end
