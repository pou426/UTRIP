c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,omega0,qt0,dvy,margin,x,ix,z,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dzm(jx),z(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension zg(jgx),dzmg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0

      omega0=1.0d-3
      qt0=1.5d0
      xrgn=2.
      zrgn=1.
      dvy=-qt0*omega0*xrgn

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=xrgn/dble(igx-margin*2)
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
c      dz,z

      dz0=zrgn/dble(jgx-margin*2)
      do j=1,jgx
         dzmg(j)=dz0
      enddo

      jzero=jgx/2+1
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
      betai=1.0d0/(1.d3)
      pr0=1.d-5
      b0=sqrt(8*pi*betai*pr0)

      do j=1,jx
      do i=1,ix
         ro(i,j) = 1.
         pr(i,j) = pr0
         bx(i,j) = 0.0
         by(i,j) = 0.0
         bz(i,j) = b0*(tanh((x(i)+0.2)/0.05)+1)/2.
     &               *(tanh((-x(i)+0.2)/0.05)+1)/2.
         vx(i,j) = 0.0
         vy(i,j) = -qt0*omega0*x(i)
         vz(i,j) = 0.0
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|

      amp=1.d-2

      mdum=-1
      do jg=1,jgx
      do ig=1,igx
         vrand=rangen(mdum)
         j=jg-jpe*(jx-2*margin)
         i=ig-ipe*(ix-2*margin)
        if (j.ge.1.and.j.le.jx .and.
     &      i.ge.1.and.i.le.ix) then
         rannum=2.d0*vrand-1.
         pr(i,j) = pr(i,j)*(1+amp*rannum)
         endif
      enddo
      enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'omega0',omega0)
      call dacputparamd(mf_params,'qt0',qt0)



      return
      end
