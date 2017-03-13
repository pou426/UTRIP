c======================================================================|
      subroutine moclagy(bycx,bycz,ro,bxm,bzm,vy,by
     &                       ,dt,dxm,dym,dzm,ix,jx,kx)
c======================================================================|
c
c NAME  moclagy
c
c PURPOSE
c    Method of characteristics (MOC) for Alfven wave 
c    for Lagrangean coordinate
c    calculate V_perp & B_perp at the cell boundary for the half step
c
c OUTPUTS
c    bycz(ix,jx,kx): [double] magnetic field at the corners
c    bycx(ix,jx,kx): [double] magnetic field at the corners
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    ro(ix,jx,kx): [double] density
c    bxm(ix,jx,kx) : [double] magnetic field
c    vy(ix,jx,kx): [double] velocity 
c    by(ix,jx,kx): [double] magnetic field
c    dxm(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c    modifide 2004-11-1 K.Takahashi
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dxm(ix),dym(jx),dzm(kx)
      dimension ro(ix,jx,kx),bxm(ix,jx,kx),bzm(ix,jx,kx)
      dimension cal(ix,jx,kx),car(ix,jx,kx)
      dimension vl(ix,jx,kx),vr(ix,jx,kx)
      dimension bl(ix,jx,kx),br(ix,jx,kx)
      dimension bycx(ix,jx,kx),bycz(ix,jx,kx)
      dimension vy(ix,jx,kx),by(ix,jx,kx)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

c-
c-  x-direction
c-

c- Characteristic velocity

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        roc=abs(0.25*(ro(i,j,k  )+ro(i+1,j,k  )
     &               +ro(i,j,k+1)+ro(i+1,j,k+1)))

        car(i,j,k) = 0.5*(bxm(i,j,k)+bxm(i,j,k+1))/sqrt(4.*pi*roc)
        cal(i,j,k) = -car(i,j,k)
      enddo
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vy,cal,car,dxm,dym,dzm,dt,ix,jx,kx,1)
      call intpvl(bl,br,by,cal,car,dxm,dym,dzm,dt,ix,jx,kx,1)

c- MOC
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        ror=0.5*(ro(i  ,j,k)+ro(i  ,j+1,k))
        rol=0.5*(ro(i+1,j,k)+ro(i+1,j+1,k))
        ror4=sqrt(abs(4.*pi*ror))
        rol4=sqrt(abs(4.*pi*rol))

        bycx(i,j,k)=(-vr(i,j,k)+vl(i,j,k)+br(i,j,k)/ror4+bl(i,j,k)/rol4)
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
        roc=abs(0.25*(ro(i,j,k  )+ro(i+1,j,k  )
     &               +ro(i,j,k+1)+ro(i+1,j,k+1)))

        car(i,j,k) = 0.5*(bzm(i,j,k)+bzm(i+1,j,k))/sqrt(4.*pi*roc)
        cal(i,j,k) = -car(i,j,k)
      enddo
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vy,cal,car,dxm,dym,dzm,dt,ix,jx,kx,3)
      call intpvl(bl,br,by,cal,car,dxm,dym,dzm,dt,ix,jx,kx,3)

c- MOC
      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        ror=0.5*(ro(i,j,k  )+ro(i,j+1,k  ))
        rol=0.5*(ro(i,j,k+1)+ro(i,j+1,k+1))
        ror4=sqrt(abs(4.*pi*ror))
        rol4=sqrt(abs(4.*pi*rol))

        bycz(i,j,k)=(-vr(i,j,k)+vl(i,j,k)+br(i,j,k)/ror4+bl(i,j,k)/rol4)
     &             /(1./ror4+1./rol4)
      enddo
      enddo
      enddo

      return
      end
