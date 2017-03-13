c======================================================================|
      subroutine model(mbnd,ro,pr,vx,vy,vz,bx,by,bz,gm,t
     &     ,gx,gxm,gy,gym,gz,gzm
     &       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm
     &     ,margin,ix,jx,kx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dx(ix),xm(ix),dxm(ix)
      dimension y(jx),dy(jx),ym(jx),dym(jx)
      dimension z(kx),dz(kx),zm(kx),dzm(kx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension gx(ix,jx,kx),gxm(ix,jx,kx)
      dimension gy(ix,jx,kx),gym(ix,jx,kx)
      dimension gz(ix,jx,kx),gzm(ix,jx,kx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=4.d0/3.d0

c----------------------------------------------------------------------|
c   for boundary condition
c----------------------------------------------------------------------|
      if (mbnd.eq.1) goto 1000
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      xmin=1.d0
      xmax=5.d0

      dx0=(xmax-xmin)/dble(max(ix-margin*2,1))
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      if (izero.lt.ix) then
        x(izero)=xmin+dxm(izero)/2.d0
      else
        izero=1
        x(izero)=0.d0
      endif
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      call grdrdy2(dx,xm,dxm,x,ix)

      ymin=0.02d0*pi
      ymax=0.98d0*pi

      dy0=(ymax-ymin)/dble(max(jx-margin*2,1))
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=margin+1
      if (jzero.lt.jx) then
        y(jzero)=ymin+dym(jzero)/2.d0
      else
        jzero=1
        y(jzero)=0.d0
      endif
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

      call grdrdy2(dy,ym,dym,y,jx)

      dz0=2.d0*pi/dble(max(kx-margin*2,1))
      do k=1,kx
         dzm(k)=dz0
      enddo

      kzero=margin+1
      if (kzero.lt.kx) then
        z(kzero)=dzm(kzero)/2.d0
      else
        kzero=1
        z(kzero)=0.d0
      endif
      do k=kzero+1,kx
         z(k) = z(k-1)+dzm(k-1)
      enddo
      do k=kzero-1,1,-1
         z(k) = z(k+1)-dzm(k)
      enddo

      call grdrdy2(dz,zm,dzm,z,kx)

c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
        gx(i,j,k) = -1.d0/x(i)**2
        gxm(i,j,k)= -1.d0/xm(i)**2
        gy(i,j,k) = 0.d0
        gym(i,j,k)= 0.d0
        gz(i,j,k) = 0.d0
        gzm(i,j,k)= 0.d0
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
1000  continue
c--- physical constant
c     g0: Gravitational Constant
c     rm_0: Mass of proton

      g0= 6.67d-8
      rm_0=1.673d-24

c--- normalization units
c     rMsun: Mass of the star
c     rn0_0: Number density at the surface
c     rr0  : Radius of the star

      rMsun=2.0d33
      rn0_0=1.d8
      rr0= 1.d11

      vv0= sqrt(g0*rMsun/rr0)
      tt0= rr0/vv0
      ro0= rm_0*rn0_0
      pr0= ro0*vv0**2
      bb0= sqrt(pr0)

c--- Parameters taken from Stone & Norman (1992)
c--- lambda may be wrong in S&N. So I changed it.
c     zetac_0  : "position" of outflow region
c     rnu_0    : entropy
c     eta_0    : scale of t
c     a0_0     : magnetic flux
c     rlambda_0: "position" of CME front
c     tstart_0 : Time of initial setup

      zetac_0=1.104d11
      rnu_0=2.42d20
      eta_0=5.24d-8
      a0_0=1.5d21
      rlambda_0=5.54d-11

      tstart_0=8.74d3

      zetac = zetac_0/rr0
      rnu = rnu_0/(vv0**2/ro0**(1/3.d0))
      eta = eta_0*tt0**2
      a0 = a0_0/(sqrt(pr0)*rr0**2)
      rlambda = rlambda_0*rr0
      tstart = tstart_0/tt0

      zlam=rlambda*zetac
      zlam2=sqrt(zlam)
      sj1=sqrt(2/pi/zlam)*(sin(zlam)/zlam-cos(zlam))
      p0= -zlam2*sj1

      if (mbnd.ne.1) t=tstart
      phi=sqrt(eta)*t
      rc=zetac*phi
      rs=phi**(7.d0/6.d0)
      d0=exp(-2.d0/3.d0/eta)

      if (mbnd.eq.1) then
        iend=margin
      else
        iend=ix
      endif

      do i=1,iend
        zeta=x(i)/phi
        zlam=rlambda*zeta
        zlam2=sqrt(zlam)
        sj0=sqrt(2/pi/zlam)*sin(zlam)
        sj1=sqrt(2/pi/zlam)*(sin(zlam)/zlam-cos(zlam))
      do j=1,jx
      do k=1,kx
        if (x(i).lt.rc) then
            pr(i,j,k)=1.d0/(4.d0*rnu**3*x(i)**4)
            ro(i,j,k)=1.d0/(rnu**3*x(i)**3)
            vx(i,j,k)=x(i)/t

            bx(i,j,k)=2.d0*a0/x(i)**2*(p0+zlam2*sj1)*cos(y(j))
            by(i,j,k)=-a0*rlambda/x(i)/phi*(zlam2*sj0-sj1/zlam2)
     &                *sin(y(j))
            bz(i,j,k)=a0*rlambda/x(i)/phi*(p0+zlam2*sj1)*sin(y(j))
            aa=(p0+zlam2*sj1)*sin(y(j))**2
            ro(i,j,k)=ro(i,j,k)
     &                +0.50d0/pi*a0**2/x(i)**3*p0*(4.d0-zlam**2)*aa
            pr(i,j,k)=pr(i,j,k)
     &                +0.25d0/pi*a0**2/x(i)**4*p0*(2.d0-zlam**2)*aa
        else
            bx(i,j,k)=0.d0
            by(i,j,k)=0.d0
            bz(i,j,k)=0.d0
          if (x(i).lt.rs) then
            pr(i,j,k)=(7.d0/6.d0)*d0*eta/phi**4
     &                *exp(2.d0/3.d0/eta/zeta**9)
            ro(i,j,k)=7.d0*d0/(phi**3*zeta**8)
     &                *exp(2.d0/3.d0/eta/zeta**9)
            vx(i,j,k)=x(i)/t
          else
            pr(i,j,k)=1.d-20
            ro(i,j,k)=d0/(x(i)**(26.d0/7.d0))
     &          *exp(2.d0/3.d0/eta/x(i)**(9.d0/7.d0))
            vx(i,j,k)=0.d0
          endif
        endif

        vy(i,j,k)=0.0d0
        vz(i,j,k)=0.0d0

      enddo
      enddo
      enddo

      if (mbnd.eq.1) return


c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'zetac_0',zetac_0)
      call dacputparamd(mf_params,'rnu_0',rnu_0)
      call dacputparamd(mf_params,'eta_0',eta_0)
      call dacputparamd(mf_params,'a0_0',a0_0)
      call dacputparamd(mf_params,'rlambda_0',rlambda_0)
      call dacputparamd(mf_params,'rMsun',rMsun)
      call dacputparamd(mf_params,'rr0',rr0)
      call dacputparamd(mf_params,'rr0_0',rr0_0)


      
      return
      end
