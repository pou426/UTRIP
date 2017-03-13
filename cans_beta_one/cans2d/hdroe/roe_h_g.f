c======================================================================|
      subroutine roe_h_g(ro,pr,vx,vy,dt,gm,gx,gy,dx,ix,dy,jx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)

      dimension fro(ix,jx),fee(ix,jx),frx(ix,jx),fry(ix,jx)
      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx)
      dimension row(ix,jx,2),prw(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension gx(ix,jx),gy(ix,jx)

c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do j=1,jx
      do i=1,ix
         rx(i,j)=ro(i,j)*vx(i,j)
         ry(i,j)=ro(i,j)*vy(i,j)
         v2=vx(i,j)**2+vy(i,j)**2
         ee(i,j)=pr(i,j)/(gm-1.0d0) +0.5d0*ro(i,j)*v2
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed half step
c     computation of 1st order flux f(i,l)
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=ro(i,j)
         eeh(i,j)=ee(i,j)
         rxh(i,j)=rx(i,j)
         ryh(i,j)=ry(i,j)
      enddo
      enddo

c     x - direction
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix-1
         row(i,j,1)=ro(i,j)
         prw(i,j,1)=pr(i,j)
         vxw(i,j,1)=vx(i,j)
         vyw(i,j,1)=vy(i,j)
         row(i,j,2)=ro(i+1,j)
         prw(i,j,2)=pr(i+1,j)
         vxw(i,j,2)=vx(i+1,j)
         vyw(i,j,2)=vy(i+1,j)
      enddo
      enddo

      call roeflux_h(fro,fee,frx,fry,gm,row,prw,vxw,vyw,ix,jx)

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
      enddo
      enddo

c     y - direction
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix
         row(i,j,1)=ro(i,j)
         prw(i,j,1)=pr(i,j)
         vxw(i,j,1)=vx(i,j)
         vyw(i,j,1)=vy(i,j)
         row(i,j,2)=ro(i,j+1)
         prw(i,j,2)=pr(i,j+1)
         vxw(i,j,2)=vx(i,j+1)
         vyw(i,j,2)=vy(i,j+1)
      enddo
      enddo

      call roeflux_h(fro,fee,fry,frx,gm,row,prw,vyw,vxw,ix,jx)

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
      enddo
      enddo

c     source
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-1
         srx=ro(i,j)*gx(i,j)
         rxh(i,j)=rxh(i,j)+0.5d0*dt*srx
         sry=ro(i,j)*gy(i,j)
         ryh(i,j)=ryh(i,j)+0.5d0*dt*sry
         see=ro(i,j)*(vx(i,j)*gx(i,j)+vy(i,j)*gy(i,j))
         eeh(i,j)=eeh(i,j)+0.5d0*dt*see
      enddo
      enddo

c     computation of basic variables on half step

      do j=2,jx-1
      do i=2,ix-1
         vxh(i,j)=rxh(i,j)/roh(i,j)
         vyh(i,j)=ryh(i,j)/roh(i,j)
         v2=vxh(i,j)**2+vyh(i,j)**2
         prh(i,j)=(gm-1.0d0)*(eeh(i,j)-0.5d0*roh(i,j)*v2)
      enddo
      enddo
c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
c     x - direction
c----------------------------------------------------------------------|
      call tvdminmod(1,roh,row,ix,jx)
      call tvdminmod(1,prh,prw,ix,jx)
      call tvdminmod(1,vxh,vxw,ix,jx)
      call tvdminmod(1,vyh,vyw,ix,jx)

      call roeflux_h(fro,fee,frx,fry,gm,row,prw,vxw,vyw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         ee(i,j)=ee(i,j)+dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
         rx(i,j)=rx(i,j)+dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ry(i,j)=ry(i,j)+dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
      enddo
      enddo


c     y - direction
c----------------------------------------------------------------------|
      call tvdminmod(2,roh,row,ix,jx)
      call tvdminmod(2,prh,prw,ix,jx)
      call tvdminmod(2,vxh,vxw,ix,jx)
      call tvdminmod(2,vyh,vyw,ix,jx)

      call roeflux_h(fro,fee,fry,frx,gm,row,prw,vyw,vxw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         ee(i,j)=ee(i,j)+dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
         rx(i,j)=rx(i,j)+dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ry(i,j)=ry(i,j)+dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
      enddo
      enddo

c     source
c----------------------------------------------------------------------|
      do j=3,jx-2
      do i=3,ix-2
         srx=roh(i,j)*gx(i,j)
         rx(i,j)=rx(i,j)+dt*srx
         sry=roh(i,j)*gy(i,j)
         ry(i,j)=ry(i,j)+dt*sry
         see=roh(i,j)*(vxh(i,j)*gx(i,j)+vyh(i,j)*gy(i,j))
         ee(i,j)=ee(i,j)+dt*see
      enddo
      enddo

c----------------------------------------------------------------------|
c     computation of basic variables on full step
      do j=3,jx-2
      do i=3,ix-2
         vx(i,j)=rx(i,j)/ro(i,j)
         vy(i,j)=ry(i,j)/ro(i,j)
         v2=vx(i,j)**2+vy(i,j)**2
         pr(i,j)=(gm-1.0d0)*(ee(i,j)-0.5d0*ro(i,j)*v2)
      enddo
      enddo
c----------------------------------------------------------------------|
      return
      end
