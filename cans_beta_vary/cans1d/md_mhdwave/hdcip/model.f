c======================================================================|
      subroutine model(ro,pr,vx,vy,by,bx,bxm,gm,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension ro(ix),pr(ix),vx(ix)
      dimension vy(ix),by(ix)
      dimension bx(ix),bxm(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=1.d0/real(ix-margin*2)

      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=ix/2
        x(izero)=-dxm(izero)/2
        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo


c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      betai=0.5d0
      b0=sqrt(8.d0*pi/gm*betai)

      wexp=2.d-2
      amp =1.d-1

      do i=1,ix
         ro(i) = 1.d0
         pr(i) = 1.d0/gm*(1.+amp*exp(-(x(i)/wexp)**2))
c        pr(i) = 1.d0/gm
         vx(i) = 0.d0
         vy(i) = 0.d0
         bx(i) = b0
         bxm(i) = b0
         by(i) = b0+amp*exp(-(x(i)/wexp)**2)
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|


      return
      end
