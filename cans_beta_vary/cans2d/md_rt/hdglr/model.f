c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,gx,gxm,gy,gym,
     &     margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)
      dimension gx(ix,jx),gy(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)
      dimension gymz(jx),roz(jx),prz(jx)

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
c     gravitation
c----------------------------------------------------------------------|
      g0= -1.0d0/gm

      do j=1,jx
        gymz(j)=g0
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
      prz(jx)=pr1
      do j=jx-1,1,-1
         roz(j)=ro0+(ro1-ro0)*(tanh(y(j)/wtr)+1)/2
         prz(j)=prz(j+1)-roz(j)*gymz(j)*dym(j)
      enddo


      do i=1,ix
      do j=1,jx
         ro(i,j)=roz(j)
         pr(i,j)=prz(j)
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
      do j=1,jx
      do i=1,ix
         rannum=2.d0*rangen(mdum)-1.
         vy(i,j) = amp*rannum
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
      call dacputparamd(mf_params,'ro0',ro0)

      
      return
      end
