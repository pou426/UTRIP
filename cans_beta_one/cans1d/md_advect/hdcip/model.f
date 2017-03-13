c======================================================================|
      subroutine model(ro,vx,vxm,margin,x,ix,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),vx(ix),vxm(ix)

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c      dxm,x

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
      mtest=1

      if (mtest.eq.1) then
      do i=1,ix
        vx(i) =  1.d0
        vxm(i) =  1.d0
       if (x(i).ge.-0.1d0 .and. x(i).le.0.0d0) ro(i)= (x(i)+0.1d0)*10.d0
       if (x(i).ge. 0.0d0 .and. x(i).le.0.1d0) ro(i)=-(x(i)-0.1d0)*10.d0
      enddo
      else if (mtest.eq.2) then
      pi = acos(-1.0d0)
      do i=1,ix
        vx(i) = 1.d0
        vxm(i) = 1.d0
        ro(i)=cos(2*pi*x(i))
      enddo
      else if (mtest.eq.3) then
      pi = acos(-1.0d0)
c     theta=pi*1.00d0
c     theta=pi*0.50d0
      theta=pi*0.25d0
      cs0=1.d0
      do i=1,ix
        vx(i) = cs0
        vxm(i) = cs0
        ro(i)=cos(theta*i)
      enddo
      endif
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparami(mf_params,'mtest',mtest)
      if (mtest.eq.3) then
        call dacputparamd(mf_params,'theta',theta)
        call dacputparamd(mf_params,'cs0',cs0)
      endif

      return
      end
