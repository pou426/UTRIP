c======================================================================|
      subroutine model(ro,pr,vx,gm,margin,x,ix
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ro(ix),pr(ix),vx(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=1.4d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c      dxm,x

      dx0=1.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin
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
      ro1=0.125d0
      pr1=0.1d0
      vx1=0.d0

      do i=1,ix
         ro(i)  = 1.d0
         vx(i)  = 0.d0
         if (x(i).lt.0.1d0) then
             pr(i)  = 1000.d0
         else
           if (x(i).lt.0.9d0) then
             pr(i)  = 0.01d0
           else
             pr(i)  = 100.d0
           endif
         endif
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ro1',ro1)
      call dacputparamd(mf_params,'pr1',pr1)
      call dacputparamd(mf_params,'vx1',vx1)

      
      return
      end
