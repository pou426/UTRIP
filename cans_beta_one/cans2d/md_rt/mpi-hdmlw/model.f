c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,gx,gxm,gy,gym,
     &     margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)
      dimension gx(ix,jx),gy(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)
      dimension gymz(jgx),roz(jgx),prz(jgx)
      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/dble(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=igx/2+1
      xg(izero)=dxmg(izero)/2.d0
      do i=izero+1,igx
         xg(i) = xg(i-1)+dxmg(i-1)
      enddo
      do i=izero-1,1,-1
         xg(i) = xg(i+1)-dxmg(i)
      enddo

      do i=1,ix
         ig=ipe*(ix-2*margin)+i
         x(i)=xg(ig)
         dxm(i)=dxmg(ig)
      enddo

c-----------------------------------------------------------------------
c      dy,y

      dy0=1.d0/dble(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo
       
      jzero=jgx/2+1
      yg(jzero)=dymg(jzero)/2.d0
      do j=jzero+1,jgx
         yg(j) = yg(j-1)+dymg(j-1)
      enddo
      do j=jzero-1,1,-1
         yg(j) = yg(j+1)-dymg(j)
      enddo

      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         y(j)=yg(jg)
         dym(j)=dymg(jg)
      enddo

c----------------------------------------------------------------------|
c     gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do jg=1,jgx
        gymz(jg)=g0
      enddo

      do j=1,jx
      do i=1,ix
         gx(i,j)=0.0
         gxm(i,j)=0.0
         gy(i,j)=g0
         gym(i,j)=g0
      enddo
      enddo

c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      ro0=1.0d0
      ro1=4.00d0
      wtr=0.01d0

      pr1=ro1/gm
      prz(jgx)=pr1
      do jg=jgx-1,1,-1
         roz(jg)=ro0+(ro1-ro0)*(tanh(yg(jg)/wtr)+1)/2
         prz(jg)=prz(jg+1)-roz(jg)*gymz(jg)*dymg(jg)
      enddo

      do i=1,ix
      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         ro(i,j)=roz(jg)
         pr(i,j)=prz(jg)
         vx(i,j)=0.d0
         vy(i,j)=0.d0
      enddo
      enddo
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      rlambda=1.d0/2.d0
      amp=0.01d0

      if (.false.) then
      mdum=-1
      do jg=1,jgx
      do ig=1,igx
        vrand=rangen(mdum)
        j=jg-jpe*(jx-2*margin)
        i=ig-ipe*(ix-2*margin)
        if (j.ge.1.and.j.le.jx .and.
     &      i.ge.1.and.i.le.ix) then
          rannum=2.d0*vrand-1.d0
          vy(i,j) = amp*rannum
        endif
      enddo
      enddo

      else

      wkx=2.0*pi/rlambda
      omegai=sqrt(g0*wkx*(ro0-ro1)/(ro0+ro1))

      do j=1,jx
      do i=1,ix
         if (y(j).le.0) then
           vv0=amp*exp(wkx*y(j))*omegai
           vx(i,j)= vv0*sin(wkx*x(i))
           vy(i,j)= vv0*cos(wkx*x(i))
           ro(i,j)=ro(i,j)*(1.d0+omegai*vy(i,j)/g0)
         else
           vv0=amp*exp(-wkx*y(j))*omegai
           vx(i,j)=-vv0*sin(wkx*x(i))
           vy(i,j)= vv0*cos(wkx*x(i))
           ro(i,j)=ro(i,j)*(1.d0+omegai*vy(i,j)/g0)
         endif
      enddo
      enddo

      endif

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ro1',ro1)

      
      return
      end
