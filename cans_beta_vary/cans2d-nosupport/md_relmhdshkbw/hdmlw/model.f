c======================================================================|
      subroutine model(ro,pr,vx,vy,bx,by,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=2.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=1./real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       dxm(1) =dxm(2)
       dxm(ix)=dxm(ix-1)

      izero=ix/2
      x(izero)=-dxm(izero)/2
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=1./real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=jx/2
      y(jzero)=-dym(jzero)/2
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      b0=sqrt(8*pi/gm)*1.d-2

      thini=60./180.*pi

      wtr=0.02
      ro0=1.
      ro1=0.125
      pr0=1.d-4
      pr1=0.1d-4
      bx0= - b0*sin(thini)
      bx1= + b0*sin(thini)
      by0= + b0*cos(thini)
      by1= - b0*cos(thini)

      do j=1,jx
      do i=1,ix
         ss=x(i)*cos(thini)+y(j)*sin(thini)
         ro(i,j) = ro0+(ro1-ro0)*(1+tanh(ss/wtr))/2
         pr(i,j) = pr0+(pr1-pr0)*(1+tanh(ss/wtr))/2
         bx(i,j) = +0.75*b0*cos(thini)+ bx0+(bx1-bx0)*(1+tanh(ss/wtr))/2
         by(i,j) = +0.75*b0*sin(thini)+ by0+(by1-by0)*(1+tanh(ss/wtr))/2
         vx(i,j) = 0.0
         vy(i,j) = 0.0
      enddo
      enddo

c     do j=1,jx
c     do i=1,ix
c        if (x(i).le.0) then
c          ro(i,j)  = 1.d0
c          pr(i,j)  = 1.d-4
c          by(i,j)  = b0*1.d-2
c        else
c          ro(i,j)  = 0.125d0
c          pr(i,j)  = 0.1d-4
c          by(i,j)  = -b0*1.d-2
c        endif
c        vx(i,j) = 0.0
c        vy(i,j) = 0.0
c        bx(i,j) = b0*0.75*1.d-2
c     enddo
c     enddo

c     do j=1,jx
c     do i=1,ix
c        if (y(j).le.0) then
c          ro(i,j)  = 1.d0
c          pr(i,j)  = 1.d-4
c          bx(i,j)  = b0*1.d-2
c        else
c          ro(i,j)  = 0.125d0
c          pr(i,j)  = 0.1d-4
c          bx(i,j)  = -b0*1.d-2
c        endif
c        vx(i,j) = 0.0
c        vy(i,j) = 0.0
c        by(i,j) = b0*0.75*1.d-2
c     enddo
c     enddo

c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)
      call dacputparamd(mf_params,'thini',thini)



      return
      end
