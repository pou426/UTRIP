c======================================================================|
      subroutine roe_m(ro,pr,vx,vy,vz,bx,by,bz,dt,gm,dx,ix,dy,jx,dz,kx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c     version 1.2 (2005/02/08 Yuji SATO)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx),dz(kx)

      dimension ro(ix,jx,kx),pr(ix,jx,kx),ee(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension rx(ix,jx,kx),ry(ix,jx,kx),rz(ix,jx,kx)

      dimension fro(ix,jx,kx),fee(ix,jx,kx)
      dimension frx(ix,jx,kx),fry(ix,jx,kx),frz(ix,jx,kx)
      dimension fbx(ix,jx,kx),fby(ix,jx,kx),fbz(ix,jx,kx)

      dimension roh(ix,jx,kx),eeh(ix,jx,kx),prh(ix,jx,kx)
      dimension rxh(ix,jx,kx),ryh(ix,jx,kx),rzh(ix,jx,kx)
      dimension vxh(ix,jx,kx),vyh(ix,jx,kx),vzh(ix,jx,kx)
      dimension bxh(ix,jx,kx),byh(ix,jx,kx),bzh(ix,jx,kx)

      dimension row(ix,jx,kx,2),prw(ix,jx,kx,2)
      dimension vxw(ix,jx,kx,2),vyw(ix,jx,kx,2),vzw(ix,jx,kx,2)
      dimension bxw(ix,jx,kx,2),byw(ix,jx,kx,2),bzw(ix,jx,kx,2)

c        k=27
c        do i=22,22
c        write(6,*) i,pr(i,27,k-1),pr(i,27,k),pr(i,27,k+1)
c        write(6,*) i,ro(i,27,k-1),ro(i,27,k),ro(i,27,k+1)
c        write(6,*) i,bx(i,27,k-1),bx(i,27,k),bx(i,27,k+1)
c        write(6,*) i,by(i,27,k-1),by(i,27,k),by(i,27,k+1)
c        write(6,*) i,bz(i,27,k-1),bz(i,27,k),bz(i,27,k+1)
c        write(6,*) i,vx(i,27,k-1),vx(i,27,k),vx(i,27,k+1)
c        write(6,*) i,vy(i,27,k-1),vy(i,27,k),vy(i,27,k+1)
c        write(6,*) i,vz(i,27,k-1),vz(i,27,k),vz(i,27,k+1)
c        enddo
         write(6,*) pr(22,27,26),ro(22,27,26),vx(22,27,26)
         write(6,*) 'hello'
c----------------------------------------------------------------------|
c     numerical parameters
      pi = acos(-1.0d0)
      pi4=4.0d0*pi
      pi4i=1.0d0/pi4
      pi8i=5.0d-1*pi4i
c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do k=1,kx
      do j=1,jx
      do i=1,ix   
         rx(i,j,k)=ro(i,j,k)*vx(i,j,k)
         ry(i,j,k)=ro(i,j,k)*vy(i,j,k)
         rz(i,j,k)=ro(i,j,k)*vz(i,j,k)
         v2=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         b2=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         ee(i,j,k)=pr(i,j,k)/(gm-1.0d0) +0.5d0*ro(i,j,k)*v2 +pi8i*b2
      enddo
      enddo
      enddo
         write(6,*) ee(22,27,26),rx(22,27,26),ry(22,27,26),rz(22,27,26)
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         roh(i,j,k)=ro(i,j,k)
         eeh(i,j,k)=ee(i,j,k)
         rxh(i,j,k)=rx(i,j,k)
         ryh(i,j,k)=ry(i,j,k)
         rzh(i,j,k)=rz(i,j,k)
         bxh(i,j,k)=bx(i,j,k)
         byh(i,j,k)=by(i,j,k)
         bzh(i,j,k)=bz(i,j,k)
      enddo
      enddo
      enddo

c     x - direction
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix-1
         row(i,j,k,1)=ro(i,j,k)
         prw(i,j,k,1)=pr(i,j,k)
         vxw(i,j,k,1)=vx(i,j,k)
         vyw(i,j,k,1)=vy(i,j,k)
         vzw(i,j,k,1)=vz(i,j,k)
         bxw(i,j,k,1)=bx(i,j,k)
         byw(i,j,k,1)=by(i,j,k)
         bzw(i,j,k,1)=bz(i,j,k)
         row(i,j,k,2)=ro(i+1,j,k)
         prw(i,j,k,2)=pr(i+1,j,k)
         vxw(i,j,k,2)=vx(i+1,j,k)
         vyw(i,j,k,2)=vy(i+1,j,k)
         vzw(i,j,k,2)=vz(i+1,j,k)
         bxw(i,j,k,2)=bx(i+1,j,k)
         byw(i,j,k,2)=by(i+1,j,k)
         bzw(i,j,k,2)=bz(i+1,j,k)
      enddo
      enddo
      enddo

      call roeflux_m(fro,fee,frx,fry,frz,fby,fbz,gm
     &               ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,ix,jx,kx)

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j,k)=roh(i,j,k)
     &             +0.5d0*dt*( (fro(i-1,j,k)-fro(i,j,k))/dx(i) )
         eeh(i,j,k)=eeh(i,j,k)
     &             +0.5d0*dt*( (fee(i-1,j,k)-fee(i,j,k))/dx(i) )
         rxh(i,j,k)=rxh(i,j,k)
     &             +0.5d0*dt*( (frx(i-1,j,k)-frx(i,j,k))/dx(i) )
         ryh(i,j,k)=ryh(i,j,k)
     &             +0.5d0*dt*( (fry(i-1,j,k)-fry(i,j,k))/dx(i) )
         rzh(i,j,k)=rzh(i,j,k)
     &             +0.5d0*dt*( (frz(i-1,j,k)-frz(i,j,k))/dx(i) )
         byh(i,j,k)=byh(i,j,k)
     &             +0.5d0*dt*( (fby(i-1,j,k)-fby(i,j,k))/dx(i) )
         bzh(i,j,k)=bzh(i,j,k)
     &             +0.5d0*dt*( (fbz(i-1,j,k)-fbz(i,j,k))/dx(i) )
      enddo
      enddo
      enddo

c     y - direction
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx-1
      do i=1,ix
         row(i,j,k,1)=ro(i,j,k)
         prw(i,j,k,1)=pr(i,j,k)
         vxw(i,j,k,1)=vx(i,j,k)
         vyw(i,j,k,1)=vy(i,j,k)
         vzw(i,j,k,1)=vz(i,j,k)
         bxw(i,j,k,1)=bx(i,j,k)
         byw(i,j,k,1)=by(i,j,k)
         bzw(i,j,k,1)=bz(i,j,k)
         row(i,j,k,2)=ro(i,j+1,k)
         prw(i,j,k,2)=pr(i,j+1,k)
         vxw(i,j,k,2)=vx(i,j+1,k)
         vyw(i,j,k,2)=vy(i,j+1,k)
         vzw(i,j,k,2)=vz(i,j+1,k)
         bxw(i,j,k,2)=bx(i,j+1,k)
         byw(i,j,k,2)=by(i,j+1,k)
         bzw(i,j,k,2)=bz(i,j+1,k)
      enddo
      enddo
      enddo

      call roeflux_m(fro,fee,fry,frz,frx,fbz,fbx,gm
     &               ,row,prw,vyw,vzw,vxw,byw,bzw,bxw,ix,jx,kx)

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j,k)=roh(i,j,k)
     &             +0.5d0*dt*( (fro(i,j-1,k)-fro(i,j,k))/dy(j) )
         eeh(i,j,k)=eeh(i,j,k)
     &             +0.5d0*dt*( (fee(i,j-1,k)-fee(i,j,k))/dy(j) )
         rxh(i,j,k)=rxh(i,j,k)
     &             +0.5d0*dt*( (frx(i,j-1,k)-frx(i,j,k))/dy(j) )
         ryh(i,j,k)=ryh(i,j,k)
     &             +0.5d0*dt*( (fry(i,j-1,k)-fry(i,j,k))/dy(j) )
         rzh(i,j,k)=rzh(i,j,k)
     &             +0.5d0*dt*( (frz(i,j-1,k)-frz(i,j,k))/dy(j) )
         bxh(i,j,k)=bxh(i,j,k)
     &             +0.5d0*dt*( (fbx(i,j-1,k)-fbx(i,j,k))/dy(j) )
         bzh(i,j,k)=bzh(i,j,k)
     &             +0.5d0*dt*( (fbz(i,j-1,k)-fbz(i,j,k))/dy(j) )
      enddo
      enddo
      enddo

c     z - direction
c----------------------------------------------------------------------|

      do k=1,kx-1
      do j=1,jx
      do i=1,ix
         row(i,j,k,1)=ro(i,j,k)
         prw(i,j,k,1)=pr(i,j,k)
         vxw(i,j,k,1)=vx(i,j,k)
         vyw(i,j,k,1)=vy(i,j,k)
         vzw(i,j,k,1)=vz(i,j,k)
         bxw(i,j,k,1)=bx(i,j,k)
         byw(i,j,k,1)=by(i,j,k)
         bzw(i,j,k,1)=bz(i,j,k)
         row(i,j,k,2)=ro(i,j,k+1)
         prw(i,j,k,2)=pr(i,j,k+1)
         vxw(i,j,k,2)=vx(i,j,k+1)
         vyw(i,j,k,2)=vy(i,j,k+1)
         vzw(i,j,k,2)=vz(i,j,k+1)
         bxw(i,j,k,2)=bx(i,j,k+1)
         byw(i,j,k,2)=by(i,j,k+1)
         bzw(i,j,k,2)=bz(i,j,k+1)
      enddo
      enddo
      enddo

      call roeflux_m(fro,fee,frz,frx,fry,fbx,fby,gm
     &               ,row,prw,vzw,vxw,vyw,bzw,bxw,byw,ix,jx,kx)

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j,k)=roh(i,j,k)
     &             +0.5d0*dt*( (fro(i,j,k-1)-fro(i,j,k))/dz(k) )
         eeh(i,j,k)=eeh(i,j,k)
     &             +0.5d0*dt*( (fee(i,j,k-1)-fee(i,j,k))/dz(k) )
         rxh(i,j,k)=rxh(i,j,k)
     &             +0.5d0*dt*( (frx(i,j,k-1)-frx(i,j,k))/dz(k) )
         ryh(i,j,k)=ryh(i,j,k)
     &             +0.5d0*dt*( (fry(i,j,k-1)-fry(i,j,k))/dz(k) )
         rzh(i,j,k)=rzh(i,j,k)
     &             +0.5d0*dt*( (frz(i,j,k-1)-frz(i,j,k))/dz(k) )
         bxh(i,j,k)=bxh(i,j,k)
     &             +0.5d0*dt*( (fbx(i,j,k-1)-fbx(i,j,k))/dz(k) )
         byh(i,j,k)=byh(i,j,k)
     &             +0.5d0*dt*( (fby(i,j,k-1)-fby(i,j,k))/dz(k) )
      enddo
      enddo
      enddo


      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
         vxh(i,j,k)=rxh(i,j,k)/roh(i,j,k)
         vyh(i,j,k)=ryh(i,j,k)/roh(i,j,k)
         vzh(i,j,k)=rzh(i,j,k)/roh(i,j,k)
         v2=vxh(i,j,k)**2+vyh(i,j,k)**2+vzh(i,j,k)**2
         b2=bxh(i,j,k)**2+byh(i,j,k)**2+bzh(i,j,k)**2
         prh(i,j,k)=(gm-1.0d0)*(eeh(i,j,k)-0.5d0*roh(i,j,k)*v2-pi8i*b2)
      enddo
      enddo
      enddo

c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
c     x - direction
c----------------------------------------------------------------------|
      call tvdminmod(1,roh,row,ix,jx,kx)
      call tvdminmod(1,prh,prw,ix,jx,kx)
      call tvdminmod(1,vxh,vxw,ix,jx,kx)
      call tvdminmod(1,vyh,vyw,ix,jx,kx)
      call tvdminmod(1,vzh,vzw,ix,jx,kx)
      call tvdminmod(1,bxh,bxw,ix,jx,kx)
      call tvdminmod(1,byh,byw,ix,jx,kx)
      call tvdminmod(1,bzh,bzw,ix,jx,kx)

      call roeflux_m(fro,fee,frx,fry,frz,fby,fbz,gm
     &               ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,ix,jx,kx)

      do k=3,kx-2
      do j=3,jx-2
      do i=3,ix-2
         ro(i,j,k)=ro(i,j,k)+dt*( (fro(i-1,j,k)-fro(i,j,k))/dx(i) )
         ee(i,j,k)=ee(i,j,k)+dt*( (fee(i-1,j,k)-fee(i,j,k))/dx(i) )
         rx(i,j,k)=rx(i,j,k)+dt*( (frx(i-1,j,k)-frx(i,j,k))/dx(i) )
         ry(i,j,k)=ry(i,j,k)+dt*( (fry(i-1,j,k)-fry(i,j,k))/dx(i) )
         rz(i,j,k)=rz(i,j,k)+dt*( (frz(i-1,j,k)-frz(i,j,k))/dx(i) )
         by(i,j,k)=by(i,j,k)+dt*( (fby(i-1,j,k)-fby(i,j,k))/dx(i) )
         bz(i,j,k)=bz(i,j,k)+dt*( (fbz(i-1,j,k)-fbz(i,j,k))/dx(i) )
      enddo
      enddo
      enddo

c     y - direction
c----------------------------------------------------------------------|
      call tvdminmod(2,roh,row,ix,jx,kx)
      call tvdminmod(2,prh,prw,ix,jx,kx)
      call tvdminmod(2,vxh,vxw,ix,jx,kx)
      call tvdminmod(2,vyh,vyw,ix,jx,kx)
      call tvdminmod(2,vzh,vzw,ix,jx,kx)
      call tvdminmod(2,bxh,bxw,ix,jx,kx)
      call tvdminmod(2,byh,byw,ix,jx,kx)
      call tvdminmod(2,bzh,bzw,ix,jx,kx)

      call roeflux_m(fro,fee,fry,frz,frx,fbz,fbx,gm
     &               ,row,prw,vyw,vzw,vxw,byw,bzw,bxw,ix,jx,kx)

      do k=3,kx-2
      do j=3,jx-2
      do i=3,ix-2
         ro(i,j,k)=ro(i,j,k)+dt*( (fro(i,j-1,k)-fro(i,j,k))/dy(j) )
         ee(i,j,k)=ee(i,j,k)+dt*( (fee(i,j-1,k)-fee(i,j,k))/dy(j) )
         rx(i,j,k)=rx(i,j,k)+dt*( (frx(i,j-1,k)-frx(i,j,k))/dy(j) )
         ry(i,j,k)=ry(i,j,k)+dt*( (fry(i,j-1,k)-fry(i,j,k))/dy(j) )
         rz(i,j,k)=rz(i,j,k)+dt*( (frz(i,j-1,k)-frz(i,j,k))/dy(j) )
         bx(i,j,k)=bx(i,j,k)+dt*( (fbx(i,j-1,k)-fbx(i,j,k))/dy(j) )
         bz(i,j,k)=bz(i,j,k)+dt*( (fbz(i,j-1,k)-fbz(i,j,k))/dy(j) )
      enddo
      enddo
      enddo

c     z - direction
c----------------------------------------------------------------------|
      call tvdminmod(3,roh,row,ix,jx,kx)
      call tvdminmod(3,prh,prw,ix,jx,kx)
      call tvdminmod(3,vxh,vxw,ix,jx,kx)
      call tvdminmod(3,vyh,vyw,ix,jx,kx)
      call tvdminmod(3,vzh,vzw,ix,jx,kx)
      call tvdminmod(3,bxh,bxw,ix,jx,kx)
      call tvdminmod(3,byh,byw,ix,jx,kx)
      call tvdminmod(3,bzh,bzw,ix,jx,kx)

      call roeflux_m(fro,fee,frz,frx,fry,fbx,fby,gm
     &               ,row,prw,vzw,vxw,vyw,bzw,bxw,byw,ix,jx,kx)

      do k=3,kx-2
      do j=3,jx-2
      do i=3,ix-2
         ro(i,j,k)=ro(i,j,k)+dt*( (fro(i,j,k-1)-fro(i,j,k))/dz(k) )
         ee(i,j,k)=ee(i,j,k)+dt*( (fee(i,j,k-1)-fee(i,j,k))/dz(k) )
         rx(i,j,k)=rx(i,j,k)+dt*( (frx(i,j,k-1)-frx(i,j,k))/dz(k) )
         ry(i,j,k)=ry(i,j,k)+dt*( (fry(i,j,k-1)-fry(i,j,k))/dz(k) )
         rz(i,j,k)=rz(i,j,k)+dt*( (frz(i,j,k-1)-frz(i,j,k))/dz(k) )
         bx(i,j,k)=bx(i,j,k)+dt*( (fbx(i,j,k-1)-fbx(i,j,k))/dz(k) )
         by(i,j,k)=by(i,j,k)+dt*( (fby(i,j,k-1)-fby(i,j,k))/dz(k) )
      enddo
      enddo
      enddo
         write(6,*) 'h0',fro(22,27,25),fro(22,27,26),fro(22,27,27)

c----------------------------------------------------------------------|
c     computation of basic variables on full step
         write(6,*) ee(22,27,26),ro(22,27,26),rx(22,27,26)
      do k=3,kx-2
      do j=3,jx-2
      do i=3,ix-2
         vx(i,j,k)=rx(i,j,k)/ro(i,j,k)
         vy(i,j,k)=ry(i,j,k)/ro(i,j,k)
         vz(i,j,k)=rz(i,j,k)/ro(i,j,k)
         v2=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
         b2=bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
         pr(i,j,k)=(gm-1.0d0)*(ee(i,j,k)-0.5d0*ro(i,j,k)*v2 -pi8i*b2)
      enddo
      enddo
      enddo
c----------------------------------------------------------------------|
         write(6,*) pr(22,27,26)
c        k=27
c        do i=22,22
c        write(6,*) i,pr(i,27,k-1),pr(i,27,k),pr(i,27,k+1)
c        write(6,*) i,ro(i,27,k-1),ro(i,27,k),ro(i,27,k+1)
c        write(6,*) i,bx(i,27,k-1),bx(i,27,k),bx(i,27,k+1)
c        write(6,*) i,by(i,27,k-1),by(i,27,k),by(i,27,k+1)
c        write(6,*) i,bz(i,27,k-1),bz(i,27,k),bz(i,27,k+1)
c        write(6,*) i,vx(i,27,k-1),vx(i,27,k),vx(i,27,k+1)
c        write(6,*) i,vy(i,27,k-1),vy(i,27,k),vy(i,27,k+1)
c        write(6,*) i,vz(i,27,k-1),vz(i,27,k),vz(i,27,k+1)
c        enddo
c        write(6,*) 'hello'

      return
      end



