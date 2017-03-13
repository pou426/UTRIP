c======================================================================|
      subroutine model(ro,vx,vy,by,vz,bz,bx,bxm,cs2,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension ro(ix),vx(ix)
      dimension vy(ix),by(ix)
      dimension vz(ix),bz(ix)
      dimension bx(ix),bxm(ix)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      cs2=1.d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=1.d0/(ix-margin*2)

      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin
        x(izero)=-dxm(izero)/2
        do i=izero+1,ix
           x(i) = x(i-1)+dxm(i-1)
        enddo
        do i=izero-1,1,-1
           x(i) = x(i+1)-dxm(i)
        enddo


c----------------------------------------------------------------------|
c     initial condition
c----------------------------------------------------------------------|
      ca0=sqrt(10.d0)
      ampaw=sqrt(0.9d0)
      rk=2.d0*pi
      bb0=ca0*sqrt(4.d0*pi)

      do i=1,ix
         phase=rk*x(i)
         ro(i)  = 1.d0
         vx(i)  = 0.0d0
         vy(i)  = ampaw*ca0*sin(phase)
         vz(i)  = ampaw*ca0*cos(phase)
         by(i)  = -ampaw*bb0*sin(phase)
         bz(i)  = -ampaw*bb0*cos(phase)
         bx(i)  = bb0
         bxm(i) = bb0
      enddo


c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'ca0',ca0)
      call dacputparamd(mf_params,'ampaw',ampaw)
      call dacputparamd(mf_params,'rk',rk)


      return
      end
