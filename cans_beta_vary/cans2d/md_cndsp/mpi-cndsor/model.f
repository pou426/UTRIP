c======================================================================|
      subroutine model(ro,pr,gm,rkap0,margin,x,ix,z,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension zg(jgx),dzmg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0
      rkap0=1.d0
      xmin=0.04d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/dble(igx-margin*2)
      do i=1,igx
         dxmg(i)=dx0
      enddo
       
      izero=margin+1
      xg(izero)=xmin+dxmg(izero)/2.d0
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
c      dz,z

      dz0=1.d0/dble(jgx-margin*2)
      do j=1,jgx
         dzmg(j)=dz0
      enddo
       
      jzero=margin+1
      zg(jzero)=dzmg(jzero)/2.d0
      do j=jzero+1,jgx
         zg(j) = zg(j-1)+dzmg(j-1)
      enddo
      do j=jzero-1,1,-1
         zg(j) = zg(j+1)-dzmg(j)
      enddo

      do j=1,jx
         jg=jpe*(jx-2*margin)+j
         z(j)=zg(jg)
         dzm(j)=dzmg(jg)
      enddo
c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|

      wexp=0.3d0

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.
         ss=sqrt(x(i)**2+z(j)**2)
         pr(i,j) = 1/gm*exp(-(ss/wexp)**2)
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      
      return
      end
