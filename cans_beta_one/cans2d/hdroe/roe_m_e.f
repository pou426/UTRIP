c======================================================================|
      subroutine roe_m_e(ro,pr,vx,vy,bx,by,az,dt,gm,et,dx,ix,dy,jx)
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

      dimension fro(ix,jx),fee(ix,jx),frx(ix,jx),fry(ix,jx)
      dimension fbx(ix,jx),fby(ix,jx)
      dimension roh(ix,jx),eeh(ix,jx),rxh(ix,jx),ryh(ix,jx)
      dimension prh(ix,jx),vxh(ix,jx),vyh(ix,jx),bxh(ix,jx),byh(ix,jx)
      dimension row(ix,jx,2),prw(ix,jx,2),vxw(ix,jx,2),vyw(ix,jx,2)
      dimension bxw(ix,jx,2),byw(ix,jx,2)

      dimension az(ix,jx)
      dimension cz(ix,jx),czh(ix,jx)
      dimension et(ix,jx)
      dimension ezh(ix,jx)
      dimension fgrid(ix,jx)



c----------------------------------------------------------------------|
c     numerical parameters
      pi=4.0d0*datan(1.0d0)
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

      call bbtocz(cz,bx,by,dx,dy,ix,jx)

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

      do j=1,jx-1
      do i=1,ix-1
        bec=-pi4i*(by(i  ,j)*et(i  ,j)*cz(i  ,j)
     &            +by(i+1,j)*et(i+1,j)*cz(i+1,j))/2
        fee(i,j)=fee(i,j)+bec
        ec=-(et(i  ,j)*cz(i  ,j)
     &      +et(i+1,j)*cz(i+1,j))/2
        fby(i,j)=fby(i,j)+ec
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
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

      do j=1,jx-1
      do i=1,ix-1
        bec=pi4i*(bx(i,j  )*et(i,j  )*cz(i,j  )
     &           +bx(i,j+1)*et(i,j+1)*cz(i,j+1))/2
        fee(i,j)=fee(i,j)+bec
        ec=(et(i,j  )*cz(i,j  )
     &     +et(i,j+1)*cz(i,j+1))/2
        fbx(i,j)=fbx(i,j)+ec
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
         roh(i,j)=roh(i,j)+0.5d0*dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         eeh(i,j)=eeh(i,j)+0.5d0*dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
         rxh(i,j)=rxh(i,j)+0.5d0*dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ryh(i,j)=ryh(i,j)+0.5d0*dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
         bxh(i,j)=bxh(i,j)+0.5d0*dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

c     computation of basic variables on half step
c----------------------------------------------------------------------|

      do j=2,jx-1
      do i=2,ix-1
         vxh(i,j)=rxh(i,j)/roh(i,j)
         vyh(i,j)=ryh(i,j)/roh(i,j)
         v2=vxh(i,j)**2+vyh(i,j)**2
         b2=bxh(i,j)**2+byh(i,j)**2
         prh(i,j)=(gm-1.0d0)*(eeh(i,j)-0.5d0*roh(i,j)*v2 -pi8i*b2)
      enddo
      enddo

      call bbtocz(czh,bxh,byh,dx,dy,ix,jx)

      do j=1,jx
      do i=1,ix
       ezh(i,j) = -vxh(i,j)*byh(i,j)+vyh(i,j)*bxh(i,j)+et(i,j)*czh(i,j)
      enddo
      enddo

c----------------------------------------------------------------------|
c     proceed full step
c     computation of 2nd order flux f(i,l)
c----------------------------------------------------------------------|
c     x - direction
c----------------------------------------------------------------------|
      do j=2,jx-1
      do i=2,ix-2
         row(i,j,1)=roh(i,j)
         prw(i,j,1)=prh(i,j)
         vxw(i,j,1)=vxh(i,j)
         vyw(i,j,1)=vyh(i,j)
         bxw(i,j,1)=bxh(i,j)
         byw(i,j,1)=byh(i,j)
         row(i,j,2)=roh(i+1,j)
         prw(i,j,2)=prh(i+1,j)
         vxw(i,j,2)=vxh(i+1,j)
         vyw(i,j,2)=vyh(i+1,j)
         bxw(i,j,2)=bxh(i+1,j)
         byw(i,j,2)=byh(i+1,j)
      enddo
      enddo

      call roeflux_m(fro,fee,frx,fry,fby,gm
     &               ,row,prw,vxw,vyw,bxw,byw,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        bec=-pi4i*(byh(i  ,j)*et(i  ,j)*czh(i  ,j)
     &            +byh(i+1,j)*et(i+1,j)*czh(i+1,j))/2
        fee(i,j)=fee(i,j)+bec
        ec=-(et(i  ,j)*czh(i  ,j)
     &      +et(i+1,j)*czh(i+1,j))/2
        fby(i,j)=fby(i,j)+ec
      enddo
      enddo

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vxh(i,j)
      enddo
      enddo
      call tvdminmod1(1,fro,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        v2=vxh(i,j)**2+vyh(i,j)**2
        ep = prh(i,j)*gm/(gm-1.)+0.5*roh(i,j)*v2
        fgrid(i,j)=ep*vxh(i,j) - byh(i,j)*ezh(i,j)*pi4i
      enddo
      enddo
      call tvdminmod1(1,fee,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vxh(i,j)**2+prh(i,j)
     &        +pi8i*(byh(i,j)**2-bxh(i,j)**2)
      enddo
      enddo
      call tvdminmod1(1,frx,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vyh(i,j)*vxh(i,j)-pi4i*byh(i,j)*bxh(i,j)
      enddo
      enddo
      call tvdminmod1(1,fry,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=-ezh(i,j)
      enddo
      enddo
      call tvdminmod1(1,fby,fgrid,ix,jx)



      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i-1,j)-fro(i,j))/dx(i) )
         ee(i,j)=ee(i,j)+dt*( (fee(i-1,j)-fee(i,j))/dx(i) )
         rx(i,j)=rx(i,j)+dt*( (frx(i-1,j)-frx(i,j))/dx(i) )
         ry(i,j)=ry(i,j)+dt*( (fry(i-1,j)-fry(i,j))/dx(i) )
         by(i,j)=by(i,j)+dt*( (fby(i-1,j)-fby(i,j))/dx(i) )
      enddo
      enddo


c     y - direction
c----------------------------------------------------------------------|
      do j=2,jx-2
      do i=2,ix-1
         row(i,j,1)=roh(i,j)
         prw(i,j,1)=prh(i,j)
         vxw(i,j,1)=vxh(i,j)
         vyw(i,j,1)=vyh(i,j)
         bxw(i,j,1)=bxh(i,j)
         byw(i,j,1)=byh(i,j)
         row(i,j,2)=roh(i,j+1)
         prw(i,j,2)=prh(i,j+1)
         vxw(i,j,2)=vxh(i,j+1)
         vyw(i,j,2)=vyh(i,j+1)
         bxw(i,j,2)=bxh(i,j+1)
         byw(i,j,2)=byh(i,j+1)
      enddo
      enddo


      call roeflux_m(fro,fee,fry,frx,fbx,gm
     &               ,row,prw,vyw,vxw,byw,bxw,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        bec=pi4i*(bxh(i,j  )*et(i,j  )*czh(i,j  )
     &           +bxh(i,j+1)*et(i,j+1)*czh(i,j+1))/2
        fee(i,j)=fee(i,j)+bec
        ec=(et(i,j  )*czh(i,j  )
     &     +et(i,j+1)*czh(i,j+1))/2
        fbx(i,j)=fbx(i,j)+ec
      enddo
      enddo

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vyh(i,j)
      enddo
      enddo
      call tvdminmod1(2,fro,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        v2=vxh(i,j)**2+vyh(i,j)**2
        ep = prh(i,j)*gm/(gm-1.)+0.5*roh(i,j)*v2
        fgrid(i,j)=ep*vyh(i,j) + bxh(i,j)*ezh(i,j)*pi4i
      enddo
      enddo
      call tvdminmod1(2,fee,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vxh(i,j)*vyh(i,j)-pi4i*bxh(i,j)*byh(i,j)
      enddo
      enddo
      call tvdminmod1(2,frx,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=roh(i,j)*vyh(i,j)**2+prh(i,j)
     &        +pi8i*(bxh(i,j)**2-byh(i,j)**2)
      enddo
      enddo
      call tvdminmod1(2,fry,fgrid,ix,jx)

      do j=2,jx-2
      do i=2,ix-2
        fgrid(i,j)=ezh(i,j)
      enddo
      enddo
      call tvdminmod1(2,fbx,fgrid,ix,jx)


      do j=3,jx-2
      do i=3,ix-2
         ro(i,j)=ro(i,j)+dt*( (fro(i,j-1)-fro(i,j))/dy(j) )
         ee(i,j)=ee(i,j)+dt*( (fee(i,j-1)-fee(i,j))/dy(j) )
         rx(i,j)=rx(i,j)+dt*( (frx(i,j-1)-frx(i,j))/dy(j) )
         ry(i,j)=ry(i,j)+dt*( (fry(i,j-1)-fry(i,j))/dy(j) )
         bx(i,j)=bx(i,j)+dt*( (fbx(i,j-1)-fbx(i,j))/dy(j) )
      enddo
      enddo

c     source term
c----------------------------------------------------------------------|

      do j=3,jx-2
      do i=3,ix-2
         ez=-vxh(i,j)*byh(i,j)+vyh(i,j)*bxh(i,j)+et(i,j)*czh(i,j)
         saz=-ez
         az(i,j)=az(i,j)+dt*saz
      enddo
      enddo

c     computation of basic variables on full step
c----------------------------------------------------------------------|
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
