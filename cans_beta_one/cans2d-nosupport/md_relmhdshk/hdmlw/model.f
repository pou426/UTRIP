c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm,margin,x,ix,y,jx
     &   ,tend,dtout,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)

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
      mdirection=3
      dtout=0.1d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      xmin=-1.0d0
      xmax= 1.0d0
      ymin=-1.0d0
      ymax= 1.0d0

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

      dy0=(ymax-ymin)/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin
      y(jzero)=ymin
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      if (mtest.eq.1) then
        tend=2.d0
        tend=1.d0
        vpara0   = 0.838d0
        vperp0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 10.00d0
        bperp0   = sqrt(4.0d0*pi)*18.28d0
        vpara1   =  0.620d0
        vperp1   = -0.442d0
        ro1   = 3.322d0
        pr1   = 55.33d0
        bperp1   = sqrt(4.0d0*pi)*14.49d0
        bpara0   = sqrt(4.d0*pi)*10.d0
      elseif (mtest.eq.2) then
        tend=2.5d0
        vpara0   = 0.870d0
        vperp0   = 0.056d0
        ro0   = 1.406d0
        pr0   = 2.015d0
        bperp0   = sqrt(4.0d0*pi)*2.500d0
        vpara1   = 0.650d0
        vperp1   = 0.284d0
        ro1   = 2.714d0
        pr1   = 5.135d0
        bperp1   = sqrt(4.0d0*pi)*4.520d0
        bpara0   = sqrt(4.d0*pi)*3.33d0
      elseif (mtest.eq.3) then
        tend=2.5d0
        vpara0   = 0.870d0
        vperp0   = 0.056d0
        ro0   = 1.406d0
        pr0   = 2.015d0
        bperp0   = sqrt(4.0d0*pi)*2.500d0
        vpara1   = 0.818d0
        vperp1   = 0.155d0
        ro1   = 1.725d0
        pr1   = 2.655d0
        bperp1   = sqrt(4.0d0*pi)*3.250d0
        bpara0   = sqrt(4.d0*pi)*3.33d0
      elseif (mtest.eq.4) then
        tend=1.0d0
        vpara0   =-0.894d0
        vperp0   = 0.000d0
        ro0   = 0.100d0
        pr0   = 1.000d0
        bperp0   = sqrt(4.0d0*pi)*0.000d0
        vpara1   =-0.180d0
        vperp1   =-0.500d0
        ro1   = 0.562d0
        pr1   =10.000d0
        bperp1   = sqrt(4.0d0*pi)*4.710d0
        bpara0   = sqrt(4.d0*pi)*2.00d0
      elseif (mtest.eq.5) then
        tend=2.0d0
        tend=1.5d0
        vpara0   =-0.409d0
        vperp0   =-0.740d0
        ro0   = 1.780d-3
        pr0   = 0.100d0
        bperp0   = sqrt(4.0d0*pi)*1.022d0
        vpara1   = 0.000d0
        vperp1   = 0.000d0
        ro1   = 0.010d0
        pr1   = 1.000d0
        bperp1   = sqrt(4.0d0*pi)*0.000d0
        bpara0   = sqrt(4.d0*pi)*1.00d0
      elseif (mtest.eq.6) then
        tend=25.0d0
        vpara0   = 0.000d0
        vperp0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d-4
        bperp0   = sqrt(4.0d0*pi)*0.010d0
        vpara1   = 0.000d0
        vperp1   = 0.000d0
        ro1   = 0.125d0
        pr1   = 1.000d-5
        bperp1   =-sqrt(4.0d0*pi)*0.010d0
        bpara0   = sqrt(4.d0*pi)*0.0075d0
      elseif (mtest.eq.7) then
        tend=1.0d0
        vpara0   = 0.000d0
        vperp0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d0
        bperp0   = sqrt(4.0d0*pi)*1.000d0
        vpara1   = 0.000d0
        vperp1   = 0.000d0
        ro1   = 0.125d0
        pr1   = 1.000d-1
        bperp1   =-sqrt(4.0d0*pi)*1.000d0
        bpara0   = sqrt(4.d0*pi)*0.75d0
      elseif (mtest.eq.8) then
        tend=1.0d0
        tend=0.6d0
        vpara0   = 0.000d0
        vperp0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 1.000d3
        bperp0   = sqrt(4.0d0*pi)*0.000d0
        vpara1   = 0.000d0
        vperp1   = 0.000d0
        ro1   = 0.100d0
        pr1   = 1.000d0
        bperp1   = sqrt(4.0d0*pi)*0.000d0
        bpara0   = sqrt(4.d0*pi)*1.00d0
      elseif (mtest.eq.9) then
        tend=1.0d0
        tend=0.6d0
        vpara0   = 0.000d0
        vperp0   = 0.000d0
        ro0   = 1.000d0
        pr0   = 30.00d0
        bperp0   = sqrt(4.0d0*pi)*20.00d0
        vpara1   = 0.000d0
        vperp1   = 0.000d0
        ro1   = 0.100d0
        pr1   = 1.000d0
        bperp1   = sqrt(4.0d0*pi)*0.000d0
        bpara0   = sqrt(4.d0*pi)*1.00d-10
      endif

      wtr=0.02
      if (mdirection.eq.1) then
      thini=0.d0
      do j=1,jx
      do i=1,ix
         ss=x(i)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         bx(i,j) = bpara0
         by(i,j) = bperp0+(bperp1-bperp0)*(1+tanh(ss/wtr))/2
         vx(i,j) = vpara0+(vpara1-vpara0)*(1+tanh(ss/wtr))/2
         vy(i,j) = vperp0+(vperp1-vperp0)*(1+tanh(ss/wtr))/2
      enddo
      enddo

      elseif (mdirection.eq.2) then
      thini=0.5d0*pi
      do j=1,jx
      do i=1,ix
         ss=y(j)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         bx(i,j) = bperp0+(bperp1-bperp0)*(1+tanh(ss/wtr))/2
         by(i,j) = bpara0
         vx(i,j) = vperp0+(vperp1-vperp0)*(1+tanh(ss/wtr))/2
         vy(i,j) = vpara0+(vpara1-vpara0)*(1+tanh(ss/wtr))/2
      enddo
      enddo

      elseif (mdirection.eq.3) then
      thini=60./180.*pi

      bx0= bpara0*cos(thini)-bperp0*sin(thini)
      by0= bpara0*sin(thini)+bperp0*cos(thini)
      vx0= vpara0*cos(thini)-vperp0*sin(thini)
      vy0= vpara0*sin(thini)+vperp0*cos(thini)
      bx1= bpara1*cos(thini)-bperp1*sin(thini)
      by1= bpara1*sin(thini)+bperp1*cos(thini)
      vx1= vpara1*cos(thini)-vperp1*sin(thini)
      vy1= vpara1*sin(thini)+vperp1*cos(thini)

      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         bx(i,j) = bx0+(bx1-bx0)*(1+tanh(ss/wtr))/2
         by(i,j) = by0+(by1-by0)*(1+tanh(ss/wtr))/2
         vx(i,j) = vx0+(vx1-vx0)*(1+tanh(ss/wtr))/2
         vy(i,j) = vy0+(vy1-vy0)*(1+tanh(ss/wtr))/2
      enddo
      enddo
      endif

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'thini',thini)
      call dacputparami(mf_params,'mtest',mtest)
      call dacputparami(mf_params,'mdirection',mdirection)



      return
      end
