c======================================================================|
      subroutine model(ro,pr,vx,vy,gm,margin,x,ix,y,jx
     &   ,rmach,xedge,thetain,ro0,pr0,vx0,vy0,mf_params)
c======================================================================|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|
      dimension x(ix),dxm(ix)
      dimension y(jx),dym(jx)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx)
      dimension vy(ix,jx)

c----------------------------------------------------------------------|
c   parameters
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      gm=1.4d0

c----------------------------------------------------------------------|
c   grid
c----------------------------------------------------------------------|
      dx0=4.d0/dble(ix-margin*2)
      do i=1,ix
         dxm(i)=dx0
      enddo

      izero=margin+1
      x(izero)=dxm(izero)/2.d0
      do i=izero+1,ix
         x(i) = x(i-1)+dxm(i-1)
      enddo
      do i=izero-1,1,-1
         x(i) = x(i+1)-dxm(i)
      enddo

      dy0=1.d0/dble(jx-margin*2)
      do j=1,jx
         dym(j)=dy0
      enddo

      jzero=margin+1
      y(jzero)=dym(jzero)/2.d0
      do j=jzero+1,jx
         y(j) = y(j-1)+dym(j-1)
      enddo
      do j=jzero-1,1,-1
         y(j) = y(j+1)-dym(j)
      enddo

c----------------------------------------------------------------------|
c     store initial condition into common area
c----------------------------------------------------------------------|

      xedge=1.d0/6.d0

c  initial medium
      ro0=gm
      pr0=1.d0
      vx0=0.d0
      vy0=0.d0

c  shock flow
      rmach=10.
      thetain=-30./180.*pi

      tanth=-tan(thetain)
      costh= cos(thetain)
      ro1=ro0 * ((gm+1)*rmach**2)/((gm-1)*rmach**2+2)
      pr1=pr0 * (2*gm*rmach**2-(gm-1))/(gm+1)
      cs0=sqrt(gm*pr0/ro0)
      dve= cs0*2*(rmach**2-1)/((gm+1)*rmach)
      vx1=vx0+dve*cos(thetain)
      vy1=vy0+dve*sin(thetain)


      do j=1,jx
      do i=1,ix
        if (x(i).le.-tan(thetain)*y(j)+xedge) then
          ro(i,j)  = ro1
          pr(i,j)  = pr1
          vx(i,j)  = vx1
          vy(i,j)  = vy1
        else
          ro(i,j)  = ro0
          pr(i,j)  = pr0
          vx(i,j)  = vx0
          vy(i,j)  = vy0
        endif
      enddo
      enddo
c----------------------------------------------------------------------|
c     write parameters to file
c----------------------------------------------------------------------|
      call dacputparamd(mf_params,'gm',gm)

      
      return
      end
