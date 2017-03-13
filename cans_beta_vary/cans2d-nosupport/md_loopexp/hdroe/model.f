c======================================================================|
      subroutine model(ro,pr,vx,vy,vz,bx,by,bz
     &    ,gm,margin,x,ix,y,jx
     &   ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx),bx(ix,jx),by(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=5./3.

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      pi2= 2.d0*acos(-1.0d0)

      dx0=1.d0/real(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo
       
      izero=margin+1
      x(izero)=dx0/2
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=4.d0/real(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo
       
      jzero=margin+1
      y(jzero)=dy0/2
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|
      betai=0.1
      bb0=sqrt(8*pi/gm*betai)
      do j=1,jx
      do i=1,ix
        ro(i,j)=1.d0
        pr(i,j)=1.d0/gm
        vx(i,j)=0.d0
        vy(i,j)=0.d0
        vz(i,j)=0.d0
        bx(i,j)= bb0*exp(-pi2*y(j))*cos(pi2*x(i))
        by(i,j)=-bb0*exp(-pi2*y(j))*sin(pi2*x(i))
        bz(i,j)=0.d0
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)


      return
      end
