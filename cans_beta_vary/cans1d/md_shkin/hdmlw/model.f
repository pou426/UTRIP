c======================================================================|
      subroutine model(ro,pr,vx,gm,margin,x,ix,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.0d0/3.0d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c      dxm,x

      dx0=10.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=ix/2
      x(izero)=-dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      ro1=0.5d0
      pr1=0.2d0

      ro2=1.0d0
      pr2=1.0d0

      ro3=0.1d0
      pr3=0.05d0

      do i=1,ix
         if (x(i).le.-1.d0) then
           ro(i)  = ro1
           pr(i)  = pr1
           vx(i)  = 0.d0
         else if (x(i).gt.1.d0) then
c        if (x(i).gt.1.d0) then
           ro(i)  = ro3
           pr(i)  = pr3
           vx(i)  = 0.d0
         else
           ro(i)  = ro2
           pr(i)  = pr2
           vx(i)  = 0.d0
         endif
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)
      call dacputparamd(mf_params,'vx0',vx0)
      call dacputparamd(mf_params,'vx1',vx1)

      
      return
      end
