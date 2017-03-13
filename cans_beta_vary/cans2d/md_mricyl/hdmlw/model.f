c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz,gm
     &     ,gx,gxm,gz,gzm,margin,x,ix,z,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),vz(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx)
      dimension gz(ix,jx),gzm(ix,jx)
      dimension gpot(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      xmax=3.d0
      xmin=0.04d0
      zmax=3.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=xmax/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=xmin+dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dz0=zmax/real(jx-margin*2)
      do j=1,jx
         dzm(j)=dz0
      enddo
       
      jzero=margin+1
      z(jzero)=dzm(jzero)/2.d0
      do j=jzero+1,jx
         z(j) = z(j-1)+dzm(j-1)
      enddo
      do j=jzero-1,1,-1
         z(j) = z(j+1)-dzm(j)
      enddo

c----------------------------------------------------------------------|
c   gravitation
c----------------------------------------------------------------------|
      sseps=0.2
      do i=1,ix
      do j=1,jx
         ss=x(i)    
         if (ss.gt.sseps) then
             gpot(i,j)=-1/ss
         else
           if (ss.gt.sseps/2) then
             gpot(i,j)=-(1/sseps-(ss-sseps)/sseps**2)
           else
             gpot(i,j)=-1.5/sseps
           endif
         endif
      enddo
      enddo

      do i=2,ix-1
      do j=2,jx-1
        gx(i,j)=-(gpot(i+1,j)-gpot(i-1,j))/(dxm(i-1)+dxm(i))
        gz(i,j)= 0
      enddo
      enddo
      call bdspnx(0,margin,gx,ix,jx)
      call bdsppx(0,margin,gz,ix,jx)
      call bdfrex(1,margin,gx,ix,jx)
      call bdfrex(1,margin,gz,ix,jx)
      call bdpery(0,margin,gx,ix,jx)
      call bdpery(0,margin,gz,ix,jx)
      call bdpery(1,margin,gx,ix,jx)
      call bdpery(1,margin,gz,ix,jx)

      do i=1,ix-1
      do j=1,jx-1
        gxm(i,j)=-(gpot(i+1,j)-gpot(i,j))/dxm(i)
        gzm(i,j)= 0
      enddo
      enddo
      call bdspnx(0,margin,gxm,ix,jx)
      call bdsppx(0,margin,gzm,ix,jx)
      call bdfrex(1,margin,gxm,ix,jx)
      call bdfrex(1,margin,gzm,ix,jx)
      call bdpery(0,margin,gxm,ix,jx)
      call bdpery(0,margin,gzm,ix,jx)
      call bdpery(1,margin,gxm,ix,jx)
      call bdpery(1,margin,gzm,ix,jx)
c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      aa=0.
      rn=3.
      eth=0.05
      emg=2.d-3
      tec0=1.d0
      roc0=1.d-3

      psi0=-1+1/2./(1-aa)+(rn+1)*eth
      pi = acos(-1.0d0)
      b0=sqrt(4.d0*pi*emg)

      do j=1,jx
      do i=1,ix
        ss=x(i)
        roc=roc0*exp(-gm/tec0*(gpot(i,j)+1.0))
        prc=roc*tec0/gm

        rod = 0
        prd = 0
        vyd = 0
        if (x(i).gt.sseps) then
          te= ( psi0+1/ss-1/2./(1-aa)*x(i)**(2*aa-2) ) / (rn+1)*gm
          if (te.gt.0) then
            rod = (te/gm/eth)**rn
            prd = te * rod /gm
            vyd = x(i)**(aa-1.0d0)
          endif
        else
        endif

        ro(i,j) = rod+roc
        pr(i,j) = prd+prc
        vy(i,j) = vyd
        vx(i,j) = 0
        vz(i,j) = 0
        by(i,j) = 0
        bx(i,j) = 0
        bz(i,j) = b0
      enddo
      enddo

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      amp=0.05d0
      rlambda=3.0d0/4.0d0

      do j=1,jx
      do i=1,ix
        vy(i,j) = vy(i,j)*(1+amp*sin(2*pi/rlambda*z(j)))
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'aa',aa)
      call dacputparamd(mf_params,'rn',rn)
      call dacputparamd(mf_params,'eth',eth)
      call dacputparamd(mf_params,'emg',emg)

      
      return
      end
