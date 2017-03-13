c======================================================================|
      subroutine moclagz(bznx,bzny,ro,bxm,bym,vz,bz,dt,dxm,dym,ix,jx)
c======================================================================|
c
c NAME  moclag
c
c PURPOSE
c    Method of characteristics (MOC) for Alfven wave 
c    for Lagrangean coordinate
c    calculate V_perp & B_perp at the cell boundary for the half step
c
c OUTPUTS
c    bxc(ix,jx),byc(ix,jx): [double] magnetic field at the corners
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix): [double] density
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
      dimension ro(ix,jx),bxm(ix,jx),bym(ix,jx)
      dimension cal(ix,jx),car(ix,jx)
      dimension vl(ix,jx),vr(ix,jx)
      dimension bl(ix,jx),br(ix,jx)
      dimension bznx(ix,jx),bzny(ix,jx)
      dimension vz(ix,jx),bz(ix,jx)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

c-
c-  x-direction
c-

c- Characteristic velocity

      do j=1,jx
      do i=1,ix-1
        ronx=abs((ro(i,j)+ro(i+1,j))/2)
        car(i,j)= bxm(i,j)/sqrt(4.*pi*ronx)
        cal(i,j)=-car(i,j)
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vz,cal,car,dxm,dym,dt,ix,jx,1)
      call intpvl(bl,br,bz,cal,car,dxm,dym,dt,ix,jx,1)

c- MOC
      do j=1,jx
      do i=1,ix-1
        ror4=sqrt(abs(4.*pi*ro(i,j)))
        rol4=sqrt(abs(4.*pi*ro(i+1,j)))
        bznx(i,j)=(-vr(i,j)+vl(i,j)+br(i,j)/ror4+bl(i,j)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 
      enddo 

c-
c-  y-direction
c-

c- Characteristic velocity

      do j=1,jx-1
      do i=1,ix
        rony=abs((ro(i,j)+ro(i,j+1))/2)
        car(i,j)= bym(i,j)/sqrt(4.*pi*rony)
        cal(i,j)=-car(i,j)
      enddo
      enddo

c- van Leer interpolation

      call intpvl(vl,vr,vz,cal,car,dxm,dym,dt,ix,jx,2)
      call intpvl(bl,br,bz,cal,car,dxm,dym,dt,ix,jx,2)

c- MOC
      do j=1,jx-1
      do i=1,ix
        ror4=sqrt(abs(4.*pi*ro(i,j)))
        rol4=sqrt(abs(4.*pi*ro(i,j+1)))
        bzny(i,j)=(-vr(i,j)+vl(i,j)+br(i,j)/ror4+bl(i,j)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 
      enddo 

      return
      end
