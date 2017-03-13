c======================================================================|
      subroutine model(ey,ez,by,bz,margin,x,ix,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension ey(ix),ez(ix),by(ix),bz(ix)

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
      amp1=1.d0
      wwav1=0.02
      amp2=0.5d0
      wwav2=0.01
      do i=1,ix
        ey(i)=amp1*exp(-(x(i)/wwav1)**2)
        ez(i)=0.d0
        by(i)=amp2*exp(-(x(i)/wwav2)**2)
        bz(i)=0.d0
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|

      
      return
      end
