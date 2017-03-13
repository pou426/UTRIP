c======================================================================|
      subroutine model(ro,vx,vy,by,bx,bxm,cs2,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension ro(ix),vx(ix)
      dimension vy(ix),by(ix)
      dimension bx(ix),bxm(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      cs2=1.d0

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
      b0=sqrt(8.d0*pi*cs2)

      do i=1,ix
         if (x(i).le.0) then
           ro(i)  = 1.d0
           by(i)  = b0
         else
           ro(i)  = 0.125d0
           by(i)  = -b0
         endif
         vx(i) = 0.0d0
         vy(i) = 0.0d0
         bx(i) = b0*0.75d0
         bxm(i) = b0*0.75d0
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|


      return
      end