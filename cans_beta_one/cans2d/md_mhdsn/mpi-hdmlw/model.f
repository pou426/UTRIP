c======================================================================|
      subroutine model(ro,pr,vx,vz,bx,bz,gm,margin,x,ix,z,jx
     &    ,mf_params,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension z(jx),dzm(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vz(ix,jx),bx(ix,jx),bz(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension zg(jgx),dzmg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=5.d0/3.d0
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
         pr(i,j)  = prism+(1.d0-prism)*exp(-(ss/wexp)**2)
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
