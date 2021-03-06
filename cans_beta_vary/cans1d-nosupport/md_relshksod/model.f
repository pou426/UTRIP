c======================================================================|
      subroutine model(ro,pr,vx,gm,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix)
      dimension dxm(ix)

      dimension ro(ix),pr(ix),vx(ix)


c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=4.d0/3.d0
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------
c      dxm,x

      dx0=100./real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
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
      ro1=0.125
      pr1=0.1
      vx0=0.
      vx1=0.

      do i=1,ix
         if (x(i).le.50.d0) then
           ro(i)  = 1.d5
           pr(i)  = 1.d0
           vx(i)  = 0.d0
         else
           ro(i)  = 0.125d5
           pr(i)  = 1.d-1
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
