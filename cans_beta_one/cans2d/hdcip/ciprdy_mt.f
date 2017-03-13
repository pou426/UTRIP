c======================================================================|
      subroutine ciprdy_mt(vxm,vym,bxm,bym
     &      ,rodx,vxdxm,vydxm,rody,vxdym,vydym
     &      ,ro,vx,vy,bx,by,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  ciprdy_h
c
c PURPOSE
c    derive temperature: derive velocity between grid points
c    derive gradient of physical variables
c        * hydrodynamics
c
c INPUTS & OUTPUTS
c    None
c
c OUTPUTS
c    te(ix,jx): [double] temperature
c    vxm(ix,jx): [double] velocity
c    vym(ix,jx): [double] velocity
c    rodx(ix,jx): [double] density gradient
c    tedx(ix,jx): [double] temperature gradient
c    vxdxm(ix,jx): [double] velocity gradient
c    vydxm(ix,jx): [double] velocity gradient
c    rody(ix,jx): [double] density gradient
c    tedy(ix,jx): [double] temperature gradient
c    vxdym(ix,jx): [double] velocity gradient
c    vydym(ix,jx): [double] velocity gradient
c
c INPUTS
c    NOTE: ??m(ix) is the variable array defined at grid bounds
c
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
c    vx(ix,jx): [double] velocity
c    gm: [double] polytropic index gamma
c    dx(ix),dxm(ix) : [double] grid spacing
c    dy(jx),dym(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-11-5 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)

      dimension ro(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension vxm(ix,jx),vym(ix,jx)
      dimension bxm(ix,jx),bym(ix,jx)
      dimension rodx(ix,jx),vxdxm(ix,jx),vydxm(ix,jx)
      dimension rody(ix,jx),vxdym(ix,jx),vydym(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix-1
        vxm(i,j)=(vx(i,j)+vx(i+1,j))/2.d0
        bxm(i,j)=(bx(i,j)+bx(i+1,j))/2.d0
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
        vym(i,j)=(vy(i,j)+vy(i,j+1))/2.d0
        bym(i,j)=(by(i,j)+by(i,j+1))/2.d0
      enddo
      enddo

      do j=1,jx
      do i=2,ix-1
         rodx(i,j)=( ro(i+1,j)- ro(i-1,j))/ dx(i)/2.d0
        vxdxm(i,j)=(vxm(i+1,j)-vxm(i-1,j))/dxm(i)/2.d0
        vydxm(i,j)=(vym(i+1,j)-vym(i-1,j))/dxm(i)/2.d0
      enddo
      enddo
      do j=2,jx-1
      do i=1,ix
         rody(i,j)=( ro(i,j+1)- ro(i,j-1))/ dy(j)/2.d0
        vxdym(i,j)=(vxm(i,j+1)-vxm(i,j-1))/dym(j)/2.d0
        vydym(i,j)=(vym(i,j+1)-vym(i,j-1))/dym(j)/2.d0
      enddo
      enddo


      return
      end
