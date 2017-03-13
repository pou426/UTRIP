c======================================================================|
      subroutine moclag_sm(bym,ro,bxm,vy,by,pr,glm,pm2,gm,dt,dxm,ix)
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
c    bym(ix): [double] magnetic field
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

      dimension dxm(ix)
      dimension bym(ix)
      dimension ro(ix),vy(ix),by(ix),bxm(ix)
      dimension cal(ix),car(ix)
      dimension vyl(ix),vyr(ix)
      dimension byl(ix),byr(ix)
      dimension pr(ix),glm(ix),pm2(ix)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)
      pi4 = 4.d0*pi

c- Characteristic velocity

      do i=1,ix-1
        rom=0.5*(ro(i)+ro(i+1))
        car(i)= bxm(i)/sqrt(4.*pi*rom)
        cal(i)=-bxm(i)/sqrt(4.*pi*rom)
      enddo

c- van Leer interpolation

      call intpvl(byl,byr,by,cal,car,dxm,dt,ix)

c- Characteristic velocity

      do i=1,ix-1
        rom=0.5*(ro(i)+ro(i+1))
        prm=0.5*(pr(i)+pr(i+1))
        rh=(rom+gm/(gm-1.d0)*prm)*glm(i)**2
        xi1=1.d0/(pi4*glm(i)**2)/(rh+pm2(i))
        car(i)=  xi1*bxm(i)/sqrt(4.*pi*rom)
        cal(i)= -xi1*bxm(i)/sqrt(4.*pi*rom)
      enddo

c- van Leer interpolation

      call intpvl(vyl,vyr,vy,cal,car,dxm,dt,ix)

c- MOC
      do i=1,ix-1
        ror4=sqrt(4.*pi*ro(i))
        rol4=sqrt(4.*pi*ro(i+1))
        bym(i)=(-vyr(i)+vyl(i)+byr(i)/ror4+byl(i)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 

      return
      end
