c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,margin,x,ix,y,jx,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=1.d0/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=ix/2+1
      x(izero)=dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=1.d0/dble(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=jx/2+1
      y(jzero)=dym(jzero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo
c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      ro0=1.0d0
      ro1=0.8d0
      u0=0.
      u1=0.5
      pr0=1.0/gm
      pr1=1.0/gm
      wtr=0.01

      betai=0.005d0
      b0=sqrt(8*pi/gm*betai)

      do j=1,jx
      do i=1,ix
         ro(i,j)=ro0+(ro1-ro0)*(1+tanh(y(j)/wtr))/2
         pr(i,j)=pr0+(pr1-pr0)*(1+tanh(y(j)/wtr))/2
         vx(i,j)=u0+(u1-u0)*(1+tanh(y(j)/wtr))/2
         vy(i,j)=0.d0
         vz(i,j)=0.d0
         bx(i,j)=b0
         by(i,j)=0.d0
         bz(i,j)=0.d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      rlambda=1.d0/2.d0
      amp=1.d-2

      if (.false.) then
      mdum=-1
      do j=1,jx
      do i=1,ix
         rannum=2.d0*rangen(mdum)-1.
         vy(i,j) = amp*rannum
      enddo
      enddo

      else

      alpha1=ro0/(ro0+ro1)
      alpha2=ro1/(ro0+ro1)

      wkx=2.0*pi/rlambda
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


      return
      end
