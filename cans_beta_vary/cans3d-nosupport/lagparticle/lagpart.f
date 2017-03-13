c======================================================================|
      subroutine lagpart(x,y,z,vx,vy,vz,ix,jx,kx,xp,yp,zp,mpx,dt)
c======================================================================|
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx),z(kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension xp(mpx),yp(mpx),zp(mpx)
c----------------------------------------------------------------------|      

      do mp=1,mpx
          ip=0
          jp=0
          kp=0
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
          do k=1,kx-1
             if(z(k).le.zp(mp).and.zp(mp).lt.z(k+1)) then
               kp=k
             endif
          enddo

        if((ip.ne.0).and.(jp.ne.0)) then
          rx0=(xp(mp)-x(ip))/(x(ip+1)-x(ip))
          rx1=1.0-rx0
          ry0=(yp(mp)-y(jp))/(y(jp+1)-y(jp))
          ry1=1.0-ry0
          rz0=(zp(mp)-z(kp))/(z(kp+1)-z(kp))
          rz1=1.0-rz0
          xp(mp)=xp(mp)+dt*
     &     (rz1*(rx1*(ry1*vx(ip  ,jp  ,kp  )+ry0*vx(ip  ,jp+1,kp  ))
     &          +rx0*(ry1*vx(ip+1,jp  ,kp  )+ry0*vx(ip+1,jp+1,kp  )))
     &     +rz0*(rx1*(ry1*vx(ip  ,jp  ,kp+1)+ry0*vx(ip  ,jp+1,kp+1))
     &          +rx0*(ry1*vx(ip+1,jp  ,kp+1)+ry0*vx(ip+1,jp+1,kp+1))))
          yp(mp)=yp(mp)+dt*
     &     (rz1*(rx1*(ry1*vy(ip  ,jp  ,kp  )+ry0*vy(ip  ,jp+1,kp  ))
     &          +rx0*(ry1*vy(ip+1,jp  ,kp  )+ry0*vy(ip+1,jp+1,kp  )))
     &     +rz0*(rx1*(ry1*vy(ip  ,jp  ,kp+1)+ry0*vy(ip  ,jp+1,kp+1))
     &          +rx0*(ry1*vy(ip+1,jp  ,kp+1)+ry0*vy(ip+1,jp+1,kp+1))))
          zp(mp)=zp(mp)+dt*
     &     (rz1*(rx1*(ry1*vz(ip  ,jp  ,kp  )+ry0*vz(ip  ,jp+1,kp  ))
     &          +rx0*(ry1*vz(ip+1,jp  ,kp  )+ry0*vz(ip+1,jp+1,kp  )))
     &     +rz0*(rx1*(ry1*vz(ip  ,jp  ,kp+1)+ry0*vz(ip  ,jp+1,kp+1))
     &          +rx0*(ry1*vz(ip+1,jp  ,kp+1)+ry0*vz(ip+1,jp+1,kp+1))))
       else
         xp(mp)=x(ix)*x(ix)
         yp(mp)=y(jx)*y(jx)
         zp(mp)=z(kx)*z(kx)
       endif
      enddo


      return
      end
