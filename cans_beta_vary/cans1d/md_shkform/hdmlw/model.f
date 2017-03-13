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
      pi = acos(-1.0d0)

      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
c      dxm,x

      dx0=1.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin
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
      amp=0.3d0
      rkx=2.d0
      do i=1,ix
        phase=2.d0*pi*rkx*x(i)
        ro(i)  = 1.d0   +amp*cos(phase)
        pr(i)  = 1.d0/gm+amp*cos(phase)
        vx(i)  = amp*cos(phase)
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|

      
      return
      end