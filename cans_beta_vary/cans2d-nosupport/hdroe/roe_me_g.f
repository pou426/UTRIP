c======================================================================|
      subroutine roe_me_g(mflg,ro,pr,vx,vy,bx,by,az,dt,gm
     &         ,fxrxi,fyrxi,fxryi,fyryi,roi
     &         ,gx,gy,dx,ix,dy,jx)
c======================================================================|
c     numerical solver of mhd equations by roe method with muscl
c     for ideal 1d simulation (2nd order)
c     version 1.1 (2001/08/24 naoya fukuda)
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension ee(ix,jx),rx(ix,jx),ry(ix,jx)
      dimension az(ix,jx)

      dimension fro(ix,jx),fee(ix,jx),frx(ix,jx),fry(ix,jx)
      dimension fbx(ix,jx),fby(ix,jx)
      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx),bxh(ix,jx),byh(ix,jx)
      dimension row(ix,jx,2),prw(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension bxw(ix,jx,2),byw(ix,jx,2)
      dimension fxrxi(ix,jx),fyrxi(ix,jx),fxryi(ix,jx),fyryi(ix,jx)
      dimension roi(ix,jx)
      dimension gx(ix,jx),gy(ix,jx)

c----------------------------------------------------------------------|
c     numerical parameters
      pi = acos(-1.0d0)
      pi4=4.0d0*pi
      pi4i=1.0d0/pi4
      pi8i=5.0d-1*pi4i
c----------------------------------------------------------------------|
c     computation of conservative variables w(i,l)
      do j=1,jx
      do i=1,ix
         rx(i,j)=ro(i,j)*vx(i,j)
         ry(i,j)=ro(i,j)*vy(i,j)
         v2=vx(i,j)**2+vy(i,j)**2
         b2=bx(i,j)**2+by(i,j)**2
         ee(i,j)=pr(i,j)/(gm-1.0d0) +0.5d0*ro(i,j)*v2 +pi8i*b2
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
         bxh(i,j)=bx(i,j)
         byh(i,j)=by(i,j)
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
         bxw(i,j,1)=bx(i,j)
         byw(i,j,1)=by(i,j)
         row(i,j,2)=ro(i+1,j)
         prw(i,j,2)=pr(i+1,j)
         vxw(i,j,2)=vx(i+1,j)
         vyw(i,j,2)=vy(i+1,j)
         bxw(i,j,2)=bx(i+1,j)
         byw(i,j,2)=by(i+1,j)
      enddo
      enddo

      call roeflux_m(fro,fee,frx,fry,fby,gm
     &               ,row,prw,vxw,vyw,bxw,byw,ix,jx)

      if (mflg.eq.1) then
        do j=1,jx
        do i=1,ix
          fxrxi(i,j)=frx(i,j)
          fxryi(i,j)=fry(i,j)
        enddo
        enddo
        goto 100
      endif

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
c        rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
c        ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( ((frx(i-1,j)-fxrxi(i-1,j))
     &                        -(frx(i,j)-fxrxi(i,j)))/dx(i) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( ((fry(i-1,j)-fxryi(i-1,j))
     &                        -(fry(i,j)-fxryi(i,j)))/dx(i) )
         byh(i,j)=byh(i,j)+0.5d0*dt*( (fby(i-1,j)-fby(i,j))/dx(i) )
      enddo
      enddo

100   continue

c     y - direction
c----------------------------------------------------------------------|

      do j=1,jx-1
      do i=1,ix
         row(i,j,1)=ro(i,j)
         prw(i,j,1)=pr(i,j)
         vxw(i,j,1)=vx(i,j)
         vyw(i,j,1)=vy(i,j)
         bxw(i,j,1)=bx(i,j)
         byw(i,j,1)=by(i,j)
         row(i,j,2)=ro(i,j+1)
         prw(i,j,2)=pr(i,j+1)
         vxw(i,j,2)=vx(i,j+1)
         vyw(i,j,2)=vy(i,j+1)
         bxw(i,j,2)=bx(i,j+1)
         byw(i,j,2)=by(i,j+1)
      enddo
      enddo

      call roeflux_m(fro,fee,fry,frx,fbx,gm
     &               ,row,prw,vyw,vxw,byw,bxw,ix,jx)

      if (mflg.eq.1) then
        do j=1,jx
        do i=1,ix
          fyrxi(i,j)=frx(i,j)
          fyryi(i,j)=fry(i,j)
        enddo
        enddo
        goto 200
      endif


      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
c        rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
c        ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( ((frx(i,j-1)-fyrxi(i,j-1)) 
     &                        -(frx(i,j  )-fyrxi(i,j  )))/dy(j) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( ((fry(i,j-1)-fyryi(i,j-1)) 
     &                        -(fry(i,j  )-fyryi(i,j  )))/dy(j) )
         bxh(i,j)=bxh(i,j)+0.5d0*dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

200   continue

c     source
c----------------------------------------------------------------------|

      if (mflg.eq.1) then
        do j=1,jx
        do i=1,ix
          roi(i,j)=ro(i,j)
        enddo
        enddo
      endif

      if (mflg.eq.1) return

      do j=2,jx-1
      do i=2,ix-1
c        srx=ro(i,j)*gx(i,j)
c        rxh(i,j)=rx(i,j)+0.5d0*dt*srx
c        sry=ro(i,j)*gy(i,j)
c        ryh(i,j)=ry(i,j)+0.5d0*dt*sry
         srx=(ro(i,j)-roi(i,j))*gx(i,j)
         rxh(i,j)=rx(i,j)+dt*srx
         sry=(ro(i,j)-roi(i,j))*gy(i,j)
         ryh(i,j)=ry(i,j)+dt*sry
         see=ro(i,j)*(vx(i,j)*gx(i,j)+vy(i,j)*gy(i,j))
         eeh(i,j)=ee(i,j)+0.5d0*dt*see
      enddo
      enddo

c     computation of basic variables on half step

      do j=2,jx-1
      do i=2,ix-1
         vxh(i,j)=rxh(i,j)/roh(i,j)
         vyh(i,j)=ryh(i,j)/roh(i,j)
         v2=vxh(i,j)**2+vyh(i,j)**2
         b2=bxh(i,j)**2+byh(i,j)**2
         prh(i,j)=(gm-1.0d0)*(eeh(i,j)-0.5d0*roh(i,j)*v2 -pi8i*b2)
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
      call tvdminmod(1,bxh,bxw,ix,jx)
      call tvdminmod(1,byh,byw,ix,jx)

      call roeflux_m(fro,fee,frx,fry,fby,gm
     &               ,row,prw,vxw,vyw,bxw,byw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         ee(i,j)=ee(i,j)+dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
         rx(i,j)=rx(i,j)+dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ry(i,j)=ry(i,j)+dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
c        rx(i,j)=rx(i,j)+dt*( ((frx(i-1,j)-fxrxi(i-1,j))
c    &                        -(frx(i,j)-fxrxi(i,j)))/dx(i) )
c        ry(i,j)=ry(i,j)+dt*( ((fry(i-1,j)-fxryi(i-1,j))
c    &                        -(fry(i,j)-fxryi(i,j)))/dx(i) )
         by(i,j)=by(i,j)+dt*( (fby(i-1,j)-fby(i,j))/dx(i) )
      enddo
      enddo


c     y - direction
c----------------------------------------------------------------------|
      call tvdminmod(2,roh,row,ix,jx)
      call tvdminmod(2,prh,prw,ix,jx)
      call tvdminmod(2,vxh,vxw,ix,jx)
      call tvdminmod(2,vyh,vyw,ix,jx)
      call tvdminmod(2,bxh,bxw,ix,jx)
      call tvdminmod(2,byh,byw,ix,jx)

      call roeflux_m(fro,fee,fry,frx,fbx,gm
     &               ,row,prw,vyw,vxw,byw,bxw,ix,jx)

      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         ee(i,j)=ee(i,j)+dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
         rx(i,j)=rx(i,j)+dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ry(i,j)=ry(i,j)+dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
c        rx(i,j)=rx(i,j)+dt*( ((frx(i,j-1)-fyrxi(i,j-1)) 
c    &                        -(frx(i,j  )-fyrxi(i,j  )))/dy(j) )
c        ry(i,j)=ry(i,j)+dt*( ((fry(i,j-1)-fyryi(i,j-1)) 
c    &                        -(fry(i,j  )-fyryi(i,j  )))/dy(j) )
         bx(i,j)=bx(i,j)+dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

c     source term
c----------------------------------------------------------------------|

      do j=3,jx-2
      do i=3,ix-2
c        srx=(roh(i,j)-rohi(i,j))*gx(i,j)
         srx=roh(i,j)*gx(i,j)
         rx(i,j)=rx(i,j)+dt*srx
c        sry=(roh(i,j)-rohi(i,j))*gy(i,j)
         sry=roh(i,j)*gy(i,j)
         ry(i,j)=ry(i,j)+dt*sry
         see=roh(i,j)*(vxh(i,j)*gx(i,j)+vyh(i,j)*gy(i,j))
         ee(i,j)=ee(i,j)+dt*see
      enddo
      enddo

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
         v2=vx(i,j)**2+vy(i,j)**2
         b2=bx(i,j)**2+by(i,j)**2
         pr(i,j)=(gm-1.0d0)*(ee(i,j)-0.5d0*ro(i,j)*v2 -pi8i*b2)
      enddo
      enddo
c----------------------------------------------------------------------|
      return
      end
