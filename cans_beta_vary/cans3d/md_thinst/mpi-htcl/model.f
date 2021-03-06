c======================================================================|
      subroutine model(ro,pr,gm,cool0,tecl0,decl0,margin,x,ix,y,jx,z,kx
     &           ,mf_params,igx,jgx,kgx,ipe,jpe,kpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx,kx),pr(ix,jx,kx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)
      dimension zg(kgx),dzmg(kgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      cool0=1.0d0
      tecl0=1.0d0
      decl0=1.d20
 
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/real(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=igx/2+1
      xg(izero)=dxmg(izero)/2.d0
      do i=izero+1,igx
         xg(i) = xg(i-1)+dxmg(i-1)
      enddo
      do i=izero-1,1,-1
         xg(i) = xg(i+1)-dxmg(i)
      enddo

      do i=1,ix
         ig=ipe*(ix-2*margin)+i
         x(i)=xg(ig)
         dxm(i)=dxmg(ig)
      enddo

c-----------------------------------------------------------------------
c      dy,y

      dy0=1.d0/real(jgx-margin*2)
      do j=1,jgx
         dymg(j)=dy0
      enddo
       
      jzero=jgx/2+1
      yg(jzero)=dymg(jzero)/2.d0
      do j=jzero+1,jgx
         yg(j) = yg(j-1)+dymg(j-1)
      enddo
      do j=jzero-1,1,-1
         yg(j) = yg(j+1)-dymg(j)
      enddo

      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         y(j)=yg(jg)
         dym(j)=dymg(jg)
      enddo

c-----------------------------------------------------------------------
c      dz,z

      dz0=1.d0/real(kgx-margin*2)
      do k=1,kgx
         dzmg(k)=dz0
      enddo

      kzero=kgx/2+1
      zg(kzero)=dzmg(kzero)/2.d0
      do k=kzero+1,kgx
         zg(k) = zg(k-1)+dzmg(k-1)
      enddo
      do k=kzero-1,1,-1
         zg(k) = zg(k+1)-dzmg(k)
      enddo

      do k=1,kx
         kg=kpe*(kx-2*margin)+k
         z(k)=zg(kg)
         dzm(k)=dzmg(kg)
      enddo

c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      te0=10.d0

      do k=1,kx
      do j=1,jx
      do i=1,ix
         ro(i,j,k)  = 1.d0
         pr(i,j,k) = te0/gm
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      
      return
      end
