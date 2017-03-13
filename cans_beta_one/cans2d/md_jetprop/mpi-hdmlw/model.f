c======================================================================|
      subroutine model(ro,pr,vx,vz,gm,margin,x,ix,z,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vz(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension zg(jgx),dzmg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5.d0/3.d0

      xmax=10.d0
      xmin=0.04d0
      zmax=30.d0
c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=xmax/dble(igx-margin*2)
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

      dz0=zmax/dble(jgx-margin*2)
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

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.d0
         vx(i,j) = 0.d0
         vz(i,j) = 0.d0
         pr(i,j) = 1.d0/gm
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)

      
      return
      end
