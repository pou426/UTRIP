c======================================================================|
      subroutine model(ro,pr,vx,vz,bx,bz,gm,margin,x,ix,z,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vz(ix,jx),bx(ix,jx),bz(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0
      xmin=0.04d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

      dx0=1.d0/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=xmin+dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dz0=1.d0/dble(jx-margin*2)
      do j=1,jx
         dzm(j)=dz0
      enddo
       
      jzero=margin+1
      z(jzero)=dzm(jzero)/2.d0
      do j=jzero+1,jx
         z(j) = z(j-1)+dzm(j-1)
      enddo
      do j=jzero-1,1,-1
         z(j) = z(j+1)-dzm(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      prism=1.e-8
      wexp=0.1d0

      betai=1.0d6
      b0=sqrt(prism*8*pi*betai)

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.d0
         vx(i,j) = 0.d0
         vz(i,j) = 0.d0
         ss=sqrt(x(i)**2+z(j)**2)
         pr(i,j) = prism+(1.d0-prism)*exp(-(ss/wexp)**2)
         bx(i,j) = 0.d0
         bz(i,j) = b0
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'wexp',wexp)
      call dacputparamd(mf_params,'prism',prism)
      call dacputparamd(mf_params,'betai',betai)

      
      return
      end
