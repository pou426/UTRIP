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
c  ID of the test
c  Here, K99 is Komisssarov (1999), DH03 is De villers & Hawley (2003)
c  mtest=1 : slow shock (K99)
c  mtest=2 : fast shock I  (DH03)
c  mtest=3 : fast shock II (DH03)
c  mtest=4 : switch-off fast rarefaction (K99)
c  mtest=5 : switch-off fast rarefaction (K99)
c  mtest=6 : non-relativistic Brio-Wu MHD shock tube (K99)
c  mtest=7 : relativistic Brio-Wu MHD shock tube (K99)
c  mtest=8 : shock tube I (K99)
c  mtest=9 : shock tube II (K99)

      mtest=1
      dtout=0.1d0

c----------------------------------------------------------------------|
c     grid

      if (mtest.eq.1.or.mtest.eq.2.or.mtest.eq.3) then
        xmin=-0.5d0
        xmax= 1.5d0
      elseif (mtest.eq.4) then
        xmin=-1.0d0
        xmax= 1.0d0
      elseif (mtest.eq.5) then
        xmin=-0.8d0
        xmax= 1.2d0
      elseif (mtest.eq.6.or.mtest.eq.7) then
        xmin=-1.0d0
        xmax= 1.0d0
      elseif (mtest.eq.8) then
        xmin=-1.0d0
        xmax= 1.5d0
      elseif (mtest.eq.9) then
        xmin=-1.2d0
        xmax= 1.2d0
      endif

        dx0=(xmax-xmin)/real(ix-margin*2)
        do i=1,ix
         dxm(i)=dx0
        enddo

        izero=margin
        x(izero)=xmin
        do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
        enddo

c----------------------------------------------------------------------|
c     store initial condition into common area

      if (mtest.eq.1) then
        tend=2.d0
        vx0   = 0.838d0
        vy0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 10.00d0
        by0   = sqrt(4.0d0*pi)*18.28d0
        vx1   =  0.620d0
        vy1   = -0.442d0
        ro1   = 3.322d0
        pr1   = 55.33d0
        by1   = sqrt(4.0d0*pi)*14.49d0
        bx0   = sqrt(4.d0*pi)*10.d0
      elseif (mtest.eq.2) then
        tend=2.5d0
        vx0   = 0.870d0
        vy0   = 0.056d0
        ro0   = 1.406d0
        pr0   = 2.015d0
        by0   = sqrt(4.0d0*pi)*2.500d0
        vx1   = 0.650d0
        vy1   = 0.284d0
        ro1   = 2.714d0
        pr1   = 5.135d0
        by1   = sqrt(4.0d0*pi)*4.520d0
        bx0   = sqrt(4.d0*pi)*3.33d0
      elseif (mtest.eq.3) then
        tend=2.5d0
        vx0   = 0.870d0
        vy0   = 0.056d0
        ro0   = 1.406d0
        pr0   = 2.015d0
        by0   = sqrt(4.0d0*pi)*2.500d0
        vx1   = 0.818d0
        vy1   = 0.155d0
        ro1   = 1.725d0
        pr1   = 2.655d0
        by1   = sqrt(4.0d0*pi)*3.250d0
        bx0   = sqrt(4.d0*pi)*3.33d0
      elseif (mtest.eq.4) then
        tend=1.0d0
        vx0   =-0.894d0
        vy0   = 0.000d0
        ro0   = 0.100d0
        pr0   = 1.000d0
        by0   = sqrt(4.0d0*pi)*0.000d0
        vx1   =-0.180d0
        vy1   =-0.500d0
        ro1   = 0.562d0
        pr1   =10.000d0
        by1   = sqrt(4.0d0*pi)*4.710d0
        bx0   = sqrt(4.d0*pi)*2.00d0
      elseif (mtest.eq.5) then
        tend=2.0d0
        vx0   =-0.409d0
        vy0   =-0.740d0
        ro0   = 1.780d-3
        pr0   = 0.100d0
        by0   = sqrt(4.0d0*pi)*1.022d0
        vx1   = 0.000d0
        vy1   = 0.000d0
        ro1   = 0.010d0
        pr1   = 1.000d0
        by1   = sqrt(4.0d0*pi)*0.000d0
        bx0   = sqrt(4.d0*pi)*1.00d0
      elseif (mtest.eq.6) then
        tend=25.0d0
        vx0   = 0.000d0
        vy0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d-4
        by0   = sqrt(4.0d0*pi)*0.010d0
        vx1   = 0.000d0
        vy1   = 0.000d0
        ro1   = 0.125d0
        pr1   = 1.000d-5
        by1   =-sqrt(4.0d0*pi)*0.010d0
        bx0   = sqrt(4.d0*pi)*0.0075d0
      elseif (mtest.eq.7) then
        tend=1.0d0
        vx0   = 0.000d0
        vy0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d0
        by0   = sqrt(4.0d0*pi)*1.000d0
        vx1   = 0.000d0
        vy1   = 0.000d0
        ro1   = 0.125d0
        pr1   = 1.000d-1
        by1   =-sqrt(4.0d0*pi)*1.000d0
        bx0   = sqrt(4.d0*pi)*0.75d0
      elseif (mtest.eq.8) then
        tend=1.0d0
        vx0   = 0.000d0
        vy0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d3
        by0   = sqrt(4.0d0*pi)*0.000d0
        vx1   = 0.000d0
        vy1   = 0.000d0
        ro1   = 0.100d0
        pr1   = 1.000d0
        by1   = sqrt(4.0d0*pi)*0.000d0
        bx0   = sqrt(4.d0*pi)*1.00d0
      elseif (mtest.eq.9) then
        tend=1.0d0
        vx0   = 0.000d0
        vy0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 30.00d0
        by0   = sqrt(4.0d0*pi)*20.00d0
        vx1   = 0.000d0
        vy1   = 0.000d0
        ro1   = 0.100d0
        pr1   = 1.000d0
        by1   = sqrt(4.0d0*pi)*0.000d0
        bx0   = sqrt(4.d0*pi)*1.00d-10
      endif

      xdc=0.00d0
      wdc=0.01d0
        do i=1,ix
         vx(i)=vx0+(vx1-vx0)*(tanh((x(i)-xdc)/wdc)+1.d0)/2.d0
         vy(i)=vy0+(vy1-vy0)*(tanh((x(i)-xdc)/wdc)+1.d0)/2.d0
         ro(i)=ro0+(ro1-ro0)*(tanh((x(i)-xdc)/wdc)+1.d0)/2.d0
         pr(i)=pr0+(pr1-pr0)*(tanh((x(i)-xdc)/wdc)+1.d0)/2.d0
         by(i)=by0+(by1-by0)*(tanh((x(i)-xdc)/wdc)+1.d0)/2.d0
         bx(i)  = bx0
         bxm(i) = bx0
        enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparami(mf_params,'mtest',mtest)


      return
      end
