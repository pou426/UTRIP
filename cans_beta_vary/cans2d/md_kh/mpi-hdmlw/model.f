c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi=      acos(-1.0d0)
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
c     store initial condition into common area
c----------------------------------------------------------------------|
      ro0=1.0d0
      ro1=0.8d0
      u0=0.d0
      u1=0.5d0
      pr0=1.d0/gm
      pr1=1.d0/gm
      wtr=0.01d0

      do j=1,jx
      do i=1,ix
         ro(i,j)=ro0+(ro1-ro0)*(1+tanh(y(j)/wtr))/2
         pr(i,j)=pr0+(pr1-pr0)*(1+tanh(y(j)/wtr))/2
         vx(i,j)=u0+(u1-u0)*(1+tanh(y(j)/wtr))/2
         vy(i,j)=0.d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      rlambda=1.d0/2.d0
      amp=1.d-2

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

      alpha1=ro0/(ro0+ro1)
      alpha2=ro1/(ro0+ro1)

      wkx=2*pi/rlambda
      omegar=wkx*(alpha1*u0+alpha2*u1)
      omegai=sqrt(wkx*wkx*alpha1*alpha2*(u0-u1)*(u0-u1))

      do j=1,jx
      do i=1,ix
       if (y(j).lt.0) then
         vre=amp*exp(wkx*y(j))*((wkx*u0-omegar)*cos(wkx*x(i))+
     &        omegai*sin(wkx*x(i)))
         vim=amp*exp(wkx*y(j))*((wkx*u0-omegar)*sin(wkx*x(i))-
     &        omegai*cos(wkx*x(i)))
         vy(i,j)=vy(i,j)+vre
         vx(i,j)=vx(i,j)-vim
         pr(i,j)=pr(i,j)
     &         -ro(i,j)/wkx*(omegai*vre-(wkx*u0-omegar)*vim)
       else
         vre=amp*exp(-wkx*y(j))*((wkx*u1-omegar)*cos(wkx*x(i))
     &        +omegai*sin(wkx*x(i)))
         vim=amp*exp(-wkx*y(j))*((wkx*u1-omegar)*sin(wkx*x(i))
     &        -omegai*cos(wkx*x(i)))
         vy(i,j)=vy(i,j)+vre
         vx(i,j)=vx(i,j)+vim
         pr(i,j)=pr(i,j)
     &        +ro(i,j)/wkx*(omegai*vre-(wkx*u1-omegar)*vim)
        endif
      enddo
      enddo

      endif

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'ro0',ro0)
      call dacputparamd(mf_params,'pr0',pr0)
      call dacputparamd(mf_params,'vx0',vx0)
      call dacputparamd(mf_params,'vx1',vx1)

      
      return
      end
