c======================================================================|
      subroutine lagpart(x,y,vx,vy,ix,jx,xp,yp,mpx,dt)
c======================================================================|
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension xp(mpx),yp(mpx)
c----------------------------------------------------------------------|      

      do mp=1,mpx
          ip=0
          jp=0
          do i=1,ix-1
             if(x(i).le.xp(mp).and.xp(mp).lt.x(i+1)) then
               ip=i
             endif
          enddo
          do j=1,jx-1
             if(y(j).le.yp(mp).and.yp(mp).lt.y(j+1)) then
               jp=j
             endif
          enddo

        if((ip.ne.0).and.(jp.ne.0)) then
          rx0=(xp(mp)-x(ip))/(x(ip+1)-x(ip))
          rx1=1.0-rx0
          ry0=(yp(mp)-y(jp))/(y(jp+1)-y(jp))
          ry1=1.0-ry0
          xp(mp)=xp(mp)+dt*
     &       (rx1*(ry1*vx(ip  ,jp  )+ry0*vx(ip  ,jp+1))
     &       +rx0*(ry1*vx(ip+1,jp  )+ry0*vx(ip+1,jp+1)))
          yp(mp)=yp(mp)+dt*
     &       (rx1*(ry1*vy(ip  ,jp  )+ry0*vy(ip  ,jp+1))
     &       +rx0*(ry1*vy(ip+1,jp  )+ry0*vy(ip+1,jp+1)))
       else
         xp(mp)=x(ix)*x(ix)
         yp(mp)=y(jx)*y(jx)
       endif
      enddo


      return
      end
