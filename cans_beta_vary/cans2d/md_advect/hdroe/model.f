c======================================================================|
      subroutine model(ro,vx,vxm,vy,vym,margin,x,y,ix,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),vx(ix,jx),vxm(ix,jx),vy(ix,jx),vym(ix,jx)

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=ix/2+1
      x(izero)=dxm(izero)/2
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=1.d0/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=jx/2+1
      y(jzero)=dym(jzero)/2
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         vx(i,j) = 1.d0
        vxm(i,j) = 1.d0
         vy(i,j) = 1.d0
        vym(i,j) = 1.d0
        s=sqrt(x(i)**2+y(j)**2)
        ro(i,j) = exp(-(s/0.1)**2)+1.d-5
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|

      
      return
      end
