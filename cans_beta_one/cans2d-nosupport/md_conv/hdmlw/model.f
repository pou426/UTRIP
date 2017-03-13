c======================================================================|
      subroutine model(ro,pr,vx,vy,gm
     &           ,gx,gxm,gy,gym,tem0,den0
     &     ,margin,x,ix,y,jx,rkap,rkapm,visc,viscm,hcs,hcsm,gasr)
c======================================================================|
c
c     symmetric boundaries at the top and bottom
c
c
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------|
      dimension dxm(ix),x(ix)
      dimension dym(jx),y(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)
      dimension gx(ix,jx),gxm(ix,jx)
      dimension gy(ix,jx),gym(ix,jx)

      dimension tem0(jx),pre0(jx),den0(jx),gym0(jx)
      dimension rkap(ix,jx),rkapm(ix,jx),visc(ix,jx),viscm(ix,jx)
      dimension hcs(ix,jx),hcsm(ix,jx)
c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)

c ratio of the specific heats
      gm=5.d0/3.d0
c polytrope index
      pltm=1.d0
c Prandtl number
      prnsg=1.d0
c aspect ratio of the domain
      arto=2.25d0
c depth parameter (density ratio = (lz+1)**pltm)
      lz=10.d0
c Rayleigh number at the top boundary
      rayle=1000.d0

c normalizations 
c gas constant
      gasr=1.d0/gm
c temperature gradient
      beta0=1.d0
c depth of the domain
      dp=1.d0
c 
      y0=1.d0/lz
c density at the top
      roup=1.d0
c tempature at the top
      tup=y0
c pressure at the top
      prup=roup*tup*gasr
c Specific heat at constant pressure      
      cp=gasr*gm/(gm-1.d0)
c gravitational acceleration
      g0 = (1.d0+pltm)*gasr*beta0
c Horizontal scale
      xmax=dp*arto




cc define from Rayleigh #
c      dp=-ymin
      bt=beta0-g0/cp
      rkap0=sqrt((roup*cp)**2*g0/tup*dp**4*bt/prnsg/rayle)
      visc0=rkap0*prnsg/cp

 


c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      marginx=5
      dx0=xmax/real(ix-marginx*2-1)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=marginx+1
      x(izero)=0.
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=dp/real(jx-margin*2-1)
      do j=1,jx
         dym(j)=dy0
      enddo

c      jzero=jx-margin+1
      j0=margin+1
      y(j0)=y0
      do j=j0+1,jx
        y(j) = y(j-1)+dym(j-1)
      enddo
      do j=j0-1,1,-1
        y(j) = y(j+1)-dym(j)
      enddo

c      write(6,*)'rkap ',rkap,'Prandl #',prnsg

c----------------------------------------------------------------------|
c     gravitation
c----------------------------------------------------------------------|

      do j=1,jx
         gym0(j)=g0
      do i=1,ix
         gx(i,j)=0.0
         gxm(i,j)=0.0
         gy(i,j)=g0
         gym(i,j)=g0
      enddo
      enddo

      

c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|

      do j=1,jx
         tem0(j)=beta0*y(j)
         den0(j)=roup*(y(j)/y0)**pltm
         pre0(j)=prup*(y(j)/y0)**(pltm+1.d0)
      enddo

      do i=1,ix
      do j=1,jx
         ro(i,j)=den0(j)
         pr(i,j)=pre0(j)
         vx(i,j)=0.d0
         vy(i,j)=0.d0
      enddo
      enddo

c----------------------------------------------------------------------|
c  thermal conductivity and viscosity
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         rkap(i,j)=rkap0
         visc(i,j)=visc0
         rkapm(i,j)=rkap0
         viscm(i,j)=visc0
      enddo
      enddo
      
c----------------------------------------------------------------------|
c  internal heating/cooling
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         hcs(i,j)=0.d0
         hcsm(i,j)=0.d0
      enddo
      enddo


c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|


      return
      end
