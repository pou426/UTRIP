c======================================================================|
      subroutine model(ro,pr,vx,gl,glm,gm,hh,hhm,gg,ggm,margin,x,ix
     &     ,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix)
      dimension dxm(ix)

      dimension ro(ix),pr(ix),vx(ix),gl(ix),glm(ix)

      dimension hh(ix,4),hhm(ix,4)
      dimension gg(ix,3),ggm(ix,3)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      gm=4./3.
      rstar=1.1
      rmax=5.1
c-----------------------------------------------------------------------
c     grid
c-----------------------------------------------------------------------
      dxt0=(log(rmax-1)-log(rstar-1))/real(ix-margin*2+1)

c-----------------------------------------------------------------------
c      dxm,x


      xt=log(rstar-1)
      x(1)=rstar
  
        do i=2,ix
           xt=xt+dxt0
           x(i) = exp(xt)+1.d0
        enddo
        do i=1,ix-1
           dxm(i) = x(i+1)-x(i)
        enddo
        dxm(ix)=dxm(ix-1)

c----------------------------------------------------------------------|
c       Schwarzschild metric
c----------------------------------------------------------------------|

      do i=1,ix
        alapse=sqrt(1.d0-1/x(i))
        hh(i,1)=alapse
        hh(i,2)=1.d0/alapse
        hh(i,3)=x(i)
        hh(i,4)=x(i)
        gg(i,1)=-1/(2.d0*alapse*x(i)**2)
        gg(i,2)=-alapse/x(i)
        gg(i,3)=-alapse/x(i)

        xm=x(i)+dxm(i)/2
        alapse=sqrt(1.d0-1/xm)
        hhm(i,1)=alapse
        hhm(i,2)=1.d0/alapse
        hhm(i,3)=xm
        hhm(i,4)=xm
        ggm(i,1)=-1/(2.d0*alapse*xm**2)
        ggm(i,2)=-alapse/xm
        ggm(i,3)=-alapse/xm
      enddo


c-----------------------------------------------------------------------|
c   initial temperature, density, pressure distributions
c-----------------------------------------------------------------------|
      do i=1,ix
        ro(i)=1.d0
        pr(i)=1.d-3
        vx(i)=0.d0
        gl(i)=1.d0
        glm(i)=1.d0
      enddo


      return
      end
