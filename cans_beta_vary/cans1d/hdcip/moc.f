c======================================================================|
      subroutine moc(vym,bym,ro,vxm,bxm,vy,by,dt,dxm,ix)
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

      dimension dxm(ix)
      dimension vym(ix),bym(ix)
      dimension ro(ix),vy(ix),by(ix),bxm(ix),vxm(ix)
      dimension cal(ix),car(ix)
      dimension vyl(ix),vyr(ix)
      dimension byl(ix),byr(ix)
c----------------------------------------------------------------------|

      pi = acos(-1.0d0)

c- Characteristic velocity

      romin=1.d-12
      do i=1,ix-1
c       roh=max(0.5*(ro(i)+ro(i+1)),romin)
        roh=0.5*(ro(i)+ro(i+1))
        car(i)= vxm(i)+bxm(i)/sqrt(4.*pi*roh)
        cal(i)= vxm(i)-bxm(i)/sqrt(4.*pi*roh)
      enddo

c- van Leer interpolation

      call intpvl(vyl,vyr,vy,cal,car,dxm,dt,ix)
      call intpvl(byl,byr,by,cal,car,dxm,dt,ix)

c- MOC
      do i=1,ix-1
c       ror4=sqrt(4.*pi*max(ro(i),romin))
c       rol4=sqrt(4.*pi*max(ro(i+1),romin))
        ror4=sqrt(4.*pi*ro(i))
        rol4=sqrt(4.*pi*ro(i+1))
        vym(i)=(vyr(i)*ror4+vyl(i)*rol4-byr(i)+byl(i))
     &          /(ror4+rol4)
        bym(i)=(-vyr(i)+vyl(i)+byr(i)/ror4+byl(i)/rol4)
     &          /(1./ror4+1./rol4)
      enddo 


      return
      end
