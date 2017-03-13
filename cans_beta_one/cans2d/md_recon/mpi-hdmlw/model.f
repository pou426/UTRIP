c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm
     &           ,et,etm,margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension et(ix,jx),etm(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      pi = acos(-1.0d0)

      xmax=5.d0
      ymax=20.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=xmax/dble(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=margin+1
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

      dy0=ymax/dble(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo
       
      jzero=margin+1
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
c   resistivity
c----------------------------------------------------------------------|
      etai=0.03d0
      retai=0.8d0

      do j=1,jx
      do i=1,ix
         rr=sqrt(x(i)**2+y(j)**2)
         if (rr.lt.retai) then
           et(i,j)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
         rr=sqrt((x(i)+dxm(i)/2)**2+(y(j)+dym(j)/2)**2)
         if (rr.lt.retai) then
           etm(i,j)=etai*(2*rr**3-3*retai*rr**2)/retai**3+etai
         endif
      enddo
      enddo
c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      betai=10.d0

      do j=1,jx
      do i=1,ix
         bx(i,j)=0.d0
         by(i,j)=sqrt(betai*pi*8/gm)*tanh(x(i)/0.5d0)
         pr(i,j)=1/gm*(1+betai)-by(i,j)**2/pi/8.d0
         ro(i,j)=pr(i,j)*gm

         vx(i,j)=0.0d0
         vy(i,j)=0.0d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'betai',betai)


      return
      end
