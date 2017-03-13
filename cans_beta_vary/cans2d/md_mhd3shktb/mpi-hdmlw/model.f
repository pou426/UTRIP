c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,margin,x,ix,y,jx,mf_params
     &           ,igx,jgx,ipe,jpe)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)

      dimension xg(igx),dxmg(igx)
      dimension yg(jgx),dymg(jgx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=2.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|

c-----------------------------------------------------------------------
c      dx,x

      dx0=1.d0/dble(igx-margin*2)
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

      dy0=1.d0/dble(jgx-margin*2)
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
c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      b0=sqrt(8*pi/gm)

      thini=60./180.*pi

      b00=0.75*b0
      wtr=0.02
      ro0=1.
      ro1=0.125
      pr0=1.
      pr1=0.1
      bz0= + b0
      bz1= - b0


      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         bx(i,j) = b00*cos(thini)
         by(i,j) = b00*sin(thini)
         bz(i,j) = bz0+(bz1-bz0)*(1+tanh(ss/wtr))/2
         vx(i,j) = 0.0
         vy(i,j) = 0.0
         vz(i,j) = 0.0
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'thini',thini)


      return
      end
