c======================================================================|
      subroutine moc(vxc,vyc,bxc,byc,ro,vxm,bxm,vym,bym
     &             ,dt,dxm,dym,ix,jx)
c======================================================================|
c
c NAME  moc
c
c PURPOSE
c    Method of characteristics (MOC) for Alfven wave
c    calculate V_perp & B_perp at the cell boundary for the half step
c
c OUTPUTS
c    vym(ix): [double] velocity 
c    bym(ix): [double] magnetic field
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
c    vxm(ix) : [double] velocity along the x-cordinate
c    bxm(ix) : [double] magnetic field
c    vy(ix): [double] velocity 
c    by(ix): [double] magnetic field
c    dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dxm(ix),dym(jx)
      dimension vxc(ix,jx),vyc(ix,jx),bxc(ix,jx),byc(ix,jx)
      dimension ro(ix,jx),vxm(ix,jx),vym(ix,jx),bxm(ix,jx),bym(ix,jx)
      dimension cal(ix,jx),car(ix,jx)
      dimension vl(ix,jx),vr(ix,jx)
      dimension bl(ix,jx),br(ix,jx)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

c-
c-  x-direction
c-

c- Characteristic velocity

      do j=1,jx-1
      do i=1,ix-1
        roh=abs(0.25*(ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1)))
        car(i,j)= (vxm(i,j)+vxm(i,j+1))/2
     &             +(bxm(i,j)+bxm(i,j+1))/2/sqrt(4.*pi*roh)
        cal(i,j)= (vxm(i,j)+vxm(i,j+1))/2
     &             -(bxm(i,j)+bxm(i,j+1))/2/sqrt(4.*pi*roh)
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vym,cal,car,dxm,dym,dt,ix,jx,1)
      call intpvl(bl,br,bym,cal,car,dxm,dym,dt,ix,jx,1)

c- MOC
      do j=1,jx-1
      do i=1,ix-1
        ror4=sqrt(abs(4.*pi*(ro(i,j)+ro(i,j+1))/2))
        rol4=sqrt(abs(4.*pi*(ro(i+1,j)+ro(i+1,j+1))/2))

        vyc(i,j)=(vr(i,j)*ror4+vl(i,j)*rol4-br(i,j)+bl(i,j))
     &          /(ror4+rol4)
        byc(i,j)=(-vr(i,j)+vl(i,j)+br(i,j)/ror4+bl(i,j)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 
      enddo 

c-
c-  y-direction
c-

c- Characteristic velocity

      do j=1,jx-1
      do i=1,ix-1
        roh=abs(0.25*(ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1)))
        car(i,j)= (vym(i,j)+vym(i+1,j))/2
     &             +(bym(i,j)+bym(i+1,j))/2/sqrt(4.*pi*roh)
        cal(i,j)= (vym(i,j)+vym(i+1,j))/2
     &             -(bym(i,j)+bym(i+1,j))/2/sqrt(4.*pi*roh)
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vxm,cal,car,dxm,dym,dt,ix,jx,2)
      call intpvl(bl,br,bxm,cal,car,dxm,dym,dt,ix,jx,2)

c- MOC
      do j=1,jx-1
      do i=1,ix-1
        ror4=sqrt(abs(4.*pi*(ro(i,j)+ro(i+1,j))/2))
        rol4=sqrt(abs(4.*pi*(ro(i,j+1)+ro(i+1,j+1))/2))

        vxc(i,j)=(vr(i,j)*ror4+vl(i,j)*rol4-br(i,j)+bl(i,j))
     &          /(ror4+rol4)
        bxc(i,j)=(-vr(i,j)+vl(i,j)+br(i,j)/ror4+bl(i,j)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 
      enddo 

      return
      end
