c======================================================================|
      subroutine moclagx(bxcy,bxcz,ro,bym,bzm,vx,bx
     &                       ,dt,dxm,dym,dzm,ix,jx,kx)
c======================================================================|
c
c NAME  moclagx
c
c PURPOSE
c    Method of characteristics (MOC) for Alfven wave 
c    for Lagrangean coordinate
c    calculate V_perp & B_perp at the cell boundary for the half step
c
c OUTPUTS
c    bxcy(ix,jx,kx): [double] magnetic field at the corners
c    bycz(ix,jx,kx): [double] magnetic field at the corners
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    ro(ix,jx,kx): [double] density
c    bym(ix,jx,kx) : [double] magnetic field
c    vz(ix,jx,kx): [double] velocity 
c    bz(ix,jx,kx): [double] magnetic field
c    dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c    modified 2004-11-1 K.Takahashi
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dxm(ix),dym(jx),dzm(kx)
      dimension ro(ix,jx,kx),bym(ix,jx,kx),bzm(ix,jx,kx)
      dimension cal(ix,jx,kx),car(ix,jx,kx)
      dimension vl(ix,jx,kx),vr(ix,jx,kx)
      dimension bl(ix,jx,kx),br(ix,jx,kx)
      dimension bxcy(ix,jx,kx),bxcz(ix,jx,kx)
      dimension vx(ix,jx,kx),bx(ix,jx,kx)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

c-
c-  y-direction
c-

c- Characteristic velocity

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        roc=abs(0.25*(ro(i,j  ,k)+ro(i,j+1,k  )
     &               +ro(i,j,k+1)+ro(i,j+1,k+1)))

        car(i,j,k) = 0.5*(bym(i,j,k)+bym(i,j,k+1))/sqrt(4.*pi*roc)
        cal(i,j,k) = -car(i,j,k)
      enddo
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vx,cal,car,dxm,dym,dzm,dt,ix,jx,kx,2)
      call intpvl(bl,br,bx,cal,car,dxm,dym,dzm,dt,ix,jx,kx,2)

c- MOC
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        ror=0.5*(ro(i,j  ,k)+ro(i+1,j  ,k))
        rol=0.5*(ro(i,j+1,k)+ro(i+1,j+1,k))
        ror4=sqrt(abs(4.*pi*ror))
        rol4=sqrt(abs(4.*pi*rol))

        bxcy(i,j,k)=(-vr(i,j,k)+vl(i,j,k)+br(i,j,k)/ror4+bl(i,j,k)/rol4)
     &             /(1./ror4+1./rol4)
      enddo 
      enddo
      enddo 

c-
c-  z-direction
c-

c- Characteristic velocity

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        roc=abs(0.25*(ro(i,j,k  )+ro(i,j+1,k  )
     &               +ro(i,j,k+1)+ro(i,j+1,k+1)))

        car(i,j,k) = 0.5*(bzm(i,j,k)+bzm(i,j+1,k))/sqrt(4.*pi*roc)
        cal(i,j,k) = -car(i,j,k)
      enddo
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vx,cal,car,dxm,dym,dzm,dt,ix,jx,kx,3)
      call intpvl(bl,br,bx,cal,car,dxm,dym,dzm,dt,ix,jx,kx,3)

c- MOC
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        ror=0.5*(ro(i,j,k  )+ro(i+1,j,k  ))
        rol=0.5*(ro(i,j,k+1)+ro(i+1,j,k+1))
        ror4=sqrt(abs(4.*pi*ror))
        rol4=sqrt(abs(4.*pi*rol))

        bxcz(i,j,k)=(-vr(i,j,k)+vl(i,j,k)+br(i,j,k)/ror4+bl(i,j,k)/rol4)
     &             /(1./ror4+1./rol4)
      enddo
      enddo
      enddo

      return
      end
