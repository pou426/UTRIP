c======================================================================|
      subroutine model(ro,pr,vx,vz,gm,margin,x,ix,z,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vz(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
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
c     initial condition
c----------------------------------------------------------------------|
      prism=1.e-8

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.d0
         vx(i,j) = 0.0
         vz(i,j) = 0.0
         pr(i,j) = prism
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'prism',prism)

      
      return
      end
