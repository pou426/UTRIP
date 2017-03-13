c======================================================================|
      subroutine model(u,margin,x,ix,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension u(ix)

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
      pi = acos(-1.0d0)
      do i=1,ix
        u(i) =  0.d0
       if (abs(x(i)).le.0.2d0) u(i)= -sin(2.d0*pi*x(i)/0.4d0)
      enddo
      else if (mtest.eq.2) then
      else if (mtest.eq.3) then
      endif
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparami(mf_params,'mtest',mtest)

      return
      end
