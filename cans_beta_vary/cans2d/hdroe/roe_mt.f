c======================================================================|
      subroutine roe_mt(ro,vx,vy,bx,by,az,dt,cs2,dx,ix,dy,jx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension rx(ix,jx),ry(ix,jx)
      dimension az(ix,jx)

      dimension fro(ix,jx),frx(ix,jx),fry(ix,jx)
      dimension fbx(ix,jx),fby(ix,jx)
      dimension roh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension vxh(ix,jx),vyh(ix,jx),bxh(ix,jx),byh(ix,jx)
      dimension row(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension bxw(ix,jx,2),byw(ix,jx,2)

c----------------------------------------------------------------------|
c     numerical parameters
      pi = acos(-1.0d0)
      pi4=4.0d0*pi
      pi8=8.0d0*pi
      pi4i=1.0d0/pi4
      pi8i=5.0d-1*pi4i
c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do j=1,jx
      do i=1,ix
         rx(i,j)=ro(i,j)*vx(i,j)
         ry(i,j)=ro(i,j)*vy(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=ro(i,j)
         rxh(i,j)=rx(i,j)
         ryh(i,j)=ry(i,j)
         bxh(i,j)=bx(i,j)
         byh(i,j)=by(i,j)
      enddo
      enddo

c     x - direction
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix-1
         row(i,j,1)=ro(i,j)
         vxw(i,j,1)=vx(i,j)
         vyw(i,j,1)=vy(i,j)
         bxw(i,j,1)=bx(i,j)
         byw(i,j,1)=by(i,j)
         row(i,j,2)=ro(i+1,j)
         vxw(i,j,2)=vx(i+1,j)
         vyw(i,j,2)=vy(i+1,j)
         bxw(i,j,2)=bx(i+1,j)
         byw(i,j,2)=by(i+1,j)
      enddo
      enddo

      call roeflux_mt(fro,frx,fry,fby,cs2
     &               ,row,vxw,vyw,bxw,byw,ix,jx)

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
         byh(i,j)=byh(i,j)+0.5d0*dt*( (fby(i-1,j)-fby(i,j))/dx(i) )
      enddo
      enddo

c     y - direction
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix
         row(i,j,1)=ro(i,j)
         vxw(i,j,1)=vx(i,j)
         vyw(i,j,1)=vy(i,j)
         bxw(i,j,1)=bx(i,j)
         byw(i,j,1)=by(i,j)
         row(i,j,2)=ro(i,j+1)
         vxw(i,j,2)=vx(i,j+1)
         vyw(i,j,2)=vy(i,j+1)
         bxw(i,j,2)=bx(i,j+1)
         byw(i,j,2)=by(i,j+1)
      enddo
      enddo

      call roeflux_mt(fro,fry,frx,fbx,cs2
     &               ,row,vyw,vxw,byw,bxw,ix,jx)


      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
         bxh(i,j)=bxh(i,j)+0.5d0*dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

c     computation of basic variables on half step

      do j=2,jx-1
      do i=2,ix-1
         vxh(i,j)=rxh(i,j)/roh(i,j)
         vyh(i,j)=ryh(i,j)/roh(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
c     x - direction
c----------------------------------------------------------------------|
      call tvdminmod(1,roh,row,ix,jx)
      call tvdminmod(1,vxh,vxw,ix,jx)
      call tvdminmod(1,vyh,vyw,ix,jx)
      call tvdminmod(1,bxh,bxw,ix,jx)
      call tvdminmod(1,byh,byw,ix,jx)

      call roeflux_mt(fro,frx,fry,fby,cs2
     &               ,row,vxw,vyw,bxw,byw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         rx(i,j)=rx(i,j)+dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ry(i,j)=ry(i,j)+dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
         by(i,j)=by(i,j)+dt*( (fby(i-1,j)-fby(i,j))/dx(i) )
      enddo
      enddo


c     y - direction
c----------------------------------------------------------------------|
      call tvdminmod(2,roh,row,ix,jx)
      call tvdminmod(2,vxh,vxw,ix,jx)
      call tvdminmod(2,vyh,vyw,ix,jx)
      call tvdminmod(2,bxh,bxw,ix,jx)
      call tvdminmod(2,byh,byw,ix,jx)

      call roeflux_mt(fro,fry,frx,fbx,cs2
     &               ,row,vyw,vxw,byw,bxw,ix,jx)


      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         rx(i,j)=rx(i,j)+dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ry(i,j)=ry(i,j)+dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
         bx(i,j)=bx(i,j)+dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

c     source term
c----------------------------------------------------------------------|

      do j=3,jx-2 
      do i=3,ix-2
         ez=-vxh(i,j)*byh(i,j)+vyh(i,j)*bxh(i,j)
         saz=-ez
         az(i,j)=az(i,j)+dt*saz
      enddo
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do j=3,jx-2
      do i=3,ix-2
         vx(i,j)=rx(i,j)/ro(i,j)
         vy(i,j)=ry(i,j)/ro(i,j)
      enddo
      enddo
c----------------------------------------------------------------------|
      return
      end
