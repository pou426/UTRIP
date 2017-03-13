c======================================================================|
      subroutine model(ro,pr,vx,vy,by
     &                 ,bx,bxm,gm,gx,gxm,rr,rrm,sc,scm,dv,margin,x,ix
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix) 

      dimension sc(ix),scm(ix),dv(ix)

      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)

      dimension bx(ix),bxm(ix)
      dimension gx(ix),gxm(ix)
      dimension rr(ix),rrm(ix)

      dimension te1d(ix),pr1d(ix),ro1d(ix)

c----------------------------------------------------------------------|
c     parameters
c----------------------------------------------------------------------|

      ztr=15.d0
      gm=5.d0/3.d0

c-----------------------------------------------------------------------
c   grid
c-----------------------------------------------------------------------

      dx0=0.04d0 
      i1=aint(2.0d0*(ztr)/dx0) 
      apx1=1.02d0
      dxmax1=50.d0*dx0 ! = 20

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

      izero=margin+1
        x(izero)=dxm(izero)/2
        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo

      xtop=(x(ix-margin)+dxm(ix-margin))

c----------------------------------------------------------------------|
c     parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

      gl0=1.0/gm
      ro0=1.0
      pr0=1.0/gm
      te0=1.0

      tpho=te0
      wtr=0.5
      tcor=200.0 
      
c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
! gravity is constant until x=400 (drops suddenly -> 0)

      zg0=400.d0
      wg0=10.d0
      do i=1,ix
         gx(i)  = -gl0*(1.d0-tanh((x(i)-zg0)/wg0))/2 ! equals -0.59999
      enddo
      
      do i=1,ix
         gxm(i) = -gl0*(1.d0-tanh((x(i)-zg0)/wg0))/2
      enddo

      call bdspnx(0,margin,gx,ix)
      call bdsmnx(0,margin-1,gxm,ix)

c----------------------------------------------------------------------|
c   curvature of magnetic field
c----------------------------------------------------------------------|
      do i=1,ix
        rr(i)=1.d0
        rrm(i)=1.d0
      enddo

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|
      ac = 300

      do i=1,ix
c            te1d(i)=tpho+0.5*(tcor-tpho)*(tanh((x(i)-ztr)/wtr)+1.)
            te1d(i)=tpho*(1+0.5*(ac-1)*(tanh((x(i)-ztr)/wtr)+1.))
      enddo


      i=izero+1
         ro1d(i) = (te1d(i-1)+0.5*gm*gxm(i-1)*dxm(i-1)/2)/
     &             (te1d(i)  -0.5*gm*gxm(i-1)*dxm(i-1)/2)*ro0

      do i=izero+2,ix
         ro1d(i) = (te1d(i-1)+0.5*gm*gxm(i-1)*dxm(i-1))/
     &             (te1d(i)  -0.5*gm*gxm(i-1)*dxm(i-1))*ro1d(i-1)
      enddo


      i=izero
         ro1d(i) = (te1d(i+1)-0.5*gm*gxm(i)*dxm(i)/2)/
     &             (te1d(i)  +0.5*gm*gxm(i)*dxm(i)/2)*ro0

      do i=izero-1,1,-1
         ro1d(i) = ro1d(i+1)
     &             *(te1d(i+1)-0.5*gm*gxm(i)*dxm(i))
     &             /(te1d(i)  +0.5*gm*gxm(i)*dxm(i))
      enddo

      do i=1,ix 
         pr1d(i) = pr0*(ro1d(i)/ro0)*(te1d(i)/te0)
      enddo

      
c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      onebeta=1.0d0
      
      do i=1,ix
         ro(i)  = ro1d(i)
         vx(i) = 0.0d0
         pr(i)  = pr1d(i)
c         bx(i)  = sqrt(onebeta*8*pi*pr(i)) ! this okay dont change
         vy(i) = 0.0d0
         by(i)  = 0.0d0
      enddo
      
      bx(1)=sqrt(onebeta*8*pi*pr(1))
      do i=2,ix
         bx(i) = bx(1)*(rr(1)/rr(i))**2
      enddo
      write(*,*),"bx", bx

      do i=1,ix-1
         bxm(i) = 0.5d0*(bx(i)+bx(i+1))
      enddo
      bxm(ix)=bxm(ix-1)

      
c----------------------------------------------------------------------|
c       cross section of flux tube
c----------------------------------------------------------------------|

      do i=1,ix
        sc(i)=1.d0/bx(i)
        scm(i)=1.d0/bxm(i)
      enddo
      do i=1,ix
        dv(i)=sc(i)
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ztr',ztr)
      call dacputparamd(mf_params,'tcor',tcor)
      

      return
      end
