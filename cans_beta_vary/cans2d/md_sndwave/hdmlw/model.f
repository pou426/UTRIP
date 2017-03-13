c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=ix/2+1
      x(izero)=dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo


      dy0=1.d0/dble(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=jx/2+1
      y(jzero)=dxm(izero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      wexp=1.d-2
      amp=1.d-1

      do j=1,jx
      do i=1,ix
         ss=sqrt(x(i)**2+y(j)**2)
         ro(i,j) = 1.d0*(1.d0+amp/gm*exp(-(ss/wexp)**2))
         pr(i,j) = 1.d0/gm*(1.d0+amp*exp(-(ss/wexp)**2))
         vx(i,j) = 0.0d0
         vy(i,j) = 0.0d0
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)

      
      return
      end
