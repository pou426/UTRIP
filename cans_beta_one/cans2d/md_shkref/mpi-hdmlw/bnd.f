c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,gm,time,x,y
     &     ,rmach,xedge,thetain,ro0,pr0,vx0,vy0,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)
c----------------------------------------------------------------------|      
      tanth=-tan(thetain)
      costh= cos(thetain)
      ro1=ro0 * ((gm+1)*rmach**2)/((gm-1)*rmach**2+2)
      pr1=pr0 * (2*gm*rmach**2-(gm-1))/(gm+1)
      cs0=sqrt(gm*pr0/ro0)
      dve= cs0*2*(rmach**2-1)/((gm+1)*rmach)
      vx1=vx0+dve*cos(thetain)
      vy1=vy0+dve*sin(thetain)
c----------------------------------------------------------------------|      
      if (ipe.eq.0) then
      call bdcnsx(0,margin,ro,ro1,ix,jx)
      call bdcnsx(0,margin,pr,pr1,ix,jx)
      call bdcnsx(0,margin,vx,vx1,ix,jx)
      call bdcnsx(0,margin,vy,vy1,ix,jx)
      endif

      if (ipe.eq.ipex-1) then
      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)
      endif

      if (jpe.eq.0) then

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)
      do j=1,margin
      do i=1,ix
        if (x(i).le.xedge) then
          ro(i,j)=ro1
          pr(i,j)=pr1
          vx(i,j)=vx1
          vy(i,j)=vy1
        endif
      enddo
      enddo

      endif

      if (jpe.eq.jpex-1) then

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vy,ix,jx)
      jbnd=jx-margin
      do j=1,margin
      xshk = tanth*y(jbnd+j)+xedge + cs0*rmach/costh*time
      do i=1,ix
        if (x(i).le.xshk) then
          ro(i,jbnd+j)=ro1
          pr(i,jbnd+j)=pr1
          vx(i,jbnd+j)=vx1
          vy(i,jbnd+j)=vy1
        endif
      enddo
      enddo

      endif


      return
      end
