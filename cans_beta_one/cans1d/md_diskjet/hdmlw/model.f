c======================================================================|
      subroutine model(ro,pr,vx,vy,by
     &   ,bx,bxm,gm,gx,gxm,rr,rrm,sc,scm,dv,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension vy(ix),by(ix)

      dimension sc(ix),scm(ix),dv(ix)
      dimension bx(ix),bxm(ix)
      dimension gx(ix),gxm(ix)
      dimension rr(ix),rrm(ix)
      dimension zz(ix),zzm(ix)
      dimension al(ix),alm(ix)
      dimension te(ix),tem(ix)
      dimension vym(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.2d0
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------
      dx0=10.d0/dble(ix-margin*2)
      apx1=1.004d0

c     dx0=1.d-3
c     apx1=1.03
c     dxmax1=1000*dx0
      i1=1

      do i=1,ix
         dxm(i)=dx0
      enddo

      if (i1.lt.ix) then
      do i=i1+1,ix
         dxm(i)=apx1*dxm(i-1)
c        if(dxm(i).gt.dxmax1) dxm(i)=dxmax1
      enddo
      endif
      call bdsmpx(0,margin-1,dxm,ix)

      izero=margin+1
      x(izero)=dxm(izero)/2.d0

        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo

c----------------------------------------------------------------------|
c       2-dimensional coordinate
c----------------------------------------------------------------------|
      ra=0.2d0
      a=1.5d0
      pi = acos(-1.0d0)

      rr0=1.d0
      zz0=0.d0
      al0=0.5d0*pi

      i=izero+1
      rr(i)=rr0+dxm(i-1)/2*cos(al0)
      zz(i)=zz0+dxm(i-1)/2*sin(al0)

      do i=izero+2,ix
        al(i-1)=atan(2.*a*(rr(i-1)-1.))
     &                   + 0.5*pi*exp(-(zz(i-1)/ra)**2)
        rr(i)=rr(i-1)+dxm(i-1)*cos(al(i-1))
        zz(i)=zz(i-1)+dxm(i-1)*sin(al(i-1))
      enddo
        al(ix)=atan(2.*a*(rr(ix)-1.))
     &                   + 0.5*pi*exp(-(zz(ix)/ra)*(zz(ix)/ra))

      i=izero
      rr(i)=rr0+dxm(i)/2*cos(al0)
      zz(i)=zz0-dxm(i)/2*sin(al0)
      do i=izero-1,1,-1    
        al(i+1)=atan(2.*a*(rr(i+1)-1.))
     &                   + 0.5*pi*exp(-(zz(i+1)/ra)**2)
        rr(i)=rr(i+1)+dxm(i+1)*cos(al(i+1))
        zz(i)=zz(i+1)-dxm(i+1)*sin(al(i+1))
      enddo              
        al(1)=atan(2.*a*(rr(1)-1.))
     &                   + 0.5*pi*exp(-(zz(1)/ra)*(zz(1)/ra))

      do i=1,ix-1
        rrm(i)=(rr(i)+rr(i+1))/2
        zzm(i)=(zz(i)+zz(i+1))/2
        alm(i)=(al(i)+al(i+1))/2
      enddo
        rrm(ix)=rr(ix)+dxm(ix)*cos(al(ix))/2
        zzm(ix)=zz(ix)+dxm(ix)*sin(al(ix))/2
        alm(ix)=al(ix)

c----------------------------------------------------------------------|
c       gravity
c----------------------------------------------------------------------|
      do i=1,ix
        dist2=rr(i)*rr(i)+zz(i)*zz(i)
        dist=sqrt(dist2)
        csth=rr(i)/dist
        sith=zz(i)/dist
        csal=cos(al(i))
        sial=sin(al(i))
        gx(i)=-(csal*csth+sial*sith)/dist2

        dist2=rrm(i)*rrm(i)+zzm(i)*zzm(i)
        dist=sqrt(dist2)
        csth=rrm(i)/dist
        sith=zzm(i)/dist
        csal=cos(alm(i))
        sial=sin(alm(i))
        gxm(i)=-(csal*csth+sial*sith)/dist2
      enddo
      call bdspnx(0,izero,gx,ix)
      call bdsmnx(0,izero-1,gxm,ix)

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|
      ztr=0.3d0
      wtr=0.03d0

      do i=1,ix
         fvl=sqrt(rr(i)*rr(i)/
     &             sqrt((rr(i)*rr(i)+zz(i)*zz(i))**3)
     &                )
         ff=-0.5*(tanh((zz(i)-ztr)/wtr)-1.)
         vy(i)=fvl*ff

         fvl=sqrt(rrm(i)*rrm(i)/
     &             sqrt((rrm(i)*rrm(i)+zzm(i)*zzm(i))**3)
     &                )
         ff=-0.5*(tanh((zzm(i)-ztr)/wtr)-1.)
         vym(i)=fvl*ff

         vx(i)=0.d0
      enddo

      tedsk=0.01d0
      rtecor=1.d1

      do i=1,ix
         te(i)=tedsk*(1.d0+0.5*(rtecor-1.d0)*(tanh((zz(i)-ztr)/wtr)+1.))
       tem(i)=tedsk*(1.d0+0.5*(rtecor-1.d0)*(tanh((zzm(i)-ztr)/wtr)+1.))
      enddo
      call bdsppx(0,izero,te,ix)
      call bdsmpx(0,izero-1,tem,ix)

      ro0 = 1.d0
      pr0 = 1.d0/gm * tem(izero)

      i=izero+1
         geff=gxm(i-1)+cos(alm(i-1))*vym(i-1)/rrm(i-1)
         hh=-tem(i-1)/geff/gm
         pr(i)=pr0*exp(-dxm(i-1)/2/hh)
         ro(i)=gm*pr(i)/te(i)

      do i=izero+2,ix
         geff=gxm(i-1)+cos(alm(i-1))*vym(i-1)/rrm(i-1)
         hh=-tem(i-1)/geff/gm
         pr(i)=pr(i-1)*exp(-dxm(i-1)/hh)
         ro(i)=gm*pr(i)/te(i)
      enddo
      call bdsppx(0,izero,ro,ix)
      call bdsppx(0,izero,pr,ix)

c----------------------------------------------------------------------|
c       magnetic field
c----------------------------------------------------------------------|
      betai=0.01d0
      bx0=sqrt(betai*pr0*8*pi)
      do i=1,ix
        bx(i)=bx0/rr(i)**2
        bxm(i)=bx0/rrm(i)**2
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

      
      call dacputparamd(mf_params,'tedsk',tedsk)
      call dacputparamd(mf_params,'rtecor',rtecor)
      call dacputparamd(mf_params,'betai',betai)
      call dacputparamd(mf_params,'ztr',ztr)
      mf_zz=91
      call dacdef1d(mf_zz,'zz.dac',6,ix)
      write(mf_zz) zz
      close(mf_zz)


      return
      end
