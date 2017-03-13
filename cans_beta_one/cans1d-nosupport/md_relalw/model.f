c======================================================================|
      subroutine model(ro,pr,vx,vy,by,bx,bxm,gm,margin,x,ix
     &     ,tend,dtout,mf_params)
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
      gm=4.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=3./real(ix-margin*2)

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
c     store initial condition into common area
c----------------------------------------------------------------------|
c  ID of the test
c  Here, based on De villers & Hawley (2003)
c  mtest=1 : ALF1
c  mtest=2 : ALF2
c  mtest=3 : ALF3
c  mtest=4 : ALF4

      mtest=4
      dtout=0.1d0

      if (mtest.eq.1) then
        beta=1.d-3
        vx0=0.0d0
        tend=0.9d0
      elseif (mtest.eq.2) then
        beta=1.d-1
        vx0=0.8d0
        tend=5.8d0
      elseif (mtest.eq.3) then
        beta=1.d-2
        vx0=0.1d0
        tend=1.1d0
      elseif (mtest.eq.4) then
        beta=3.15d-1
        vx0=0.2d0
        tend=5.2d0
      endif

      a0=1.d-3
      a0=1.d-5

      ro0=1.d0
      pr0=1.d-2
      bx0=sqrt(sqrt(2.d0)*4.d0*pi*pr0/beta)

      do i=1,ix
         ro(i) = ro0
         pr(i) = pr0
         vx(i) = vx0
         bx(i) = bx0
         bxm(i) = bx0
         if (x(i).le.-0.5d0.or.x(i).gt.0.5d0) then
           by(i) = 0.d0
           vy(i) = 0.d0
         else
           if (x(i).ge.-0.5d0.and.x(i).lt.0.d0) then
             by(i) = 0.d0
             vy(i) = a0
           else
             by(i) = 0.d0
             vy(i) = -a0
           endif
         endif
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|


      return
      end
