c======================================================================|
      subroutine ciprdy_m(te,vxm,vym,vzm,bxm,bym,bzm
     &                      ,rodx,tedx,vxdxm,vydxm,vzdxm
     &                      ,rody,tedy,vxdym,vydym,vzdym
     &                      ,rodz,tedz,vxdzm,vydzm,vzdzm
     &                      ,ro,pr,vx,vy,vz,bx,by,bz
     &                      ,gm,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  ciprdy_m
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
c    te(ix,jx,kx): [double] temperature
c    vxm(ix,jx,kx): [double] velocity
c    vym(ix,jx,kx): [double] velocity
c    vzm(ix,jx,kx): [double] velocity
c    rodx(ix,jx,kx): [double] density gradient
c    tedx(ix,jx,kx): [double] temperature gradient
c    vxdxm(ix,jx,kx): [double] velocity gradient
c    vydxm(ix,jx,kx): [double] velocity gradient
c    vzdxm(ix,jx,kx): [double] velocity gradient
c    rody(ix,jx,kx): [double] density gradient
c    tedy(ix,jx,kx): [double] temperature gradient
c    vxdym(ix,jx,kx): [double] velocity gradient
c    vydym(ix,jx,kx): [double] velocity gradient
c    vzdym(ix,jx,kx): [double] velocity gradient
c    rodz(ix,jx,kx): [double] density gradient
c    tedz(ix,jx,kx): [double] temperature gradient
c    vxdzm(ix,jx,kx): [double] velocity gradient
c    vydzm(ix,jx,kx): [double] velocity gradient
c    vzdzm(ix,jx,kx): [double] velocity gradient
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    ro(ix,jx,kx): [double] density
c    pr(ix,jx,kx): [double] pressure
c    vx(ix,jx,kx): [double] velocity
c    gm: [double] polytropic index gamma
c    dx(ix),dxm(ix) : [double] grid spacing
c    dy(jx),dym(jx) : [double] grid spacing
c    dz(kx),dzm(kx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)

      dimension ro(ix,jx,kx),pr(ix,jx,kx)
      dimension vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension te(ix,jx,kx),vxm(ix,jx,kx),vym(ix,jx,kx),vzm(ix,jx,kx)
      dimension bxm(ix,jx,kx),bym(ix,jx,kx),bzm(ix,jx,kx)
      dimension rodx(ix,jx,kx),rody(ix,jx,kx),rodz(ix,jx,kx)
      dimension tedx(ix,jx,kx),tedy(ix,jx,kx),tedz(ix,jx,kx)
      dimension vxdxm(ix,jx,kx),vxdym(ix,jx,kx),vxdzm(ix,jx,kx)
      dimension vydxm(ix,jx,kx),vydym(ix,jx,kx),vydzm(ix,jx,kx)
      dimension vzdxm(ix,jx,kx),vzdym(ix,jx,kx),vzdzm(ix,jx,kx)

c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
         te(i,j,k)=gm*pr(i,j,k)/ro(i,j,k)
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=1,ix-1
        vxm(i,j,k)=(vx(i,j,k)+vx(i+1,j,k))/2
        bxm(i,j,k)=(bx(i,j,k)+bx(i+1,j,k))/2
      enddo
      enddo
      enddo
      do k=1,kx
      do j=1,jx-1
      do i=1,ix
        vym(i,j,k)=(vy(i,j,k)+vy(i,j+1,k))/2
        bym(i,j,k)=(by(i,j,k)+by(i,j+1,k))/2
      enddo
      enddo
      enddo
      do k=1,kx-1
      do j=1,jx
      do i=1,ix
        vzm(i,j,k)=(vz(i,j,k)+vz(i,j,k+1))/2
        bzm(i,j,k)=(bz(i,j,k)+bz(i,j,k+1))/2
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx
      do i=2,ix-1
         rodx(i,j,k)=( ro(i+1,j,k)- ro(i-1,j,k))/ dx(i)/2
         tedx(i,j,k)=( te(i+1,j,k)- te(i-1,j,k))/ dx(i)/2
        vxdxm(i,j,k)=(vxm(i+1,j,k)-vxm(i-1,j,k))/dxm(i)/2
        vydxm(i,j,k)=(vym(i+1,j,k)-vym(i-1,j,k))/dxm(i)/2
        vzdxm(i,j,k)=(vzm(i+1,j,k)-vzm(i-1,j,k))/dxm(i)/2
      enddo
      enddo
      enddo
      do k=1,kx
      do j=2,jx-1
      do i=1,ix
         rody(i,j,k)=( ro(i,j+1,k)- ro(i,j-1,k))/ dy(j)/2
         tedy(i,j,k)=( te(i,j+1,k)- te(i,j-1,k))/ dy(j)/2
        vxdym(i,j,k)=(vxm(i,j+1,k)-vxm(i,j-1,k))/dym(j)/2
        vydym(i,j,k)=(vym(i,j+1,k)-vym(i,j-1,k))/dym(j)/2
        vzdym(i,j,k)=(vzm(i,j+1,k)-vzm(i,j-1,k))/dym(j)/2
      enddo
      enddo
      enddo
      do k=2,kx-1
      do j=1,jx
      do i=1,ix
         rodz(i,j,k)=( ro(i,j,k+1)- ro(i,j,k-1))/ dz(k)/2
         tedz(i,j,k)=( te(i,j,k+1)- te(i,j,k-1))/ dz(k)/2
        vxdzm(i,j,k)=(vxm(i,j,k+1)-vxm(i,j,k-1))/dzm(k)/2
        vydzm(i,j,k)=(vym(i,j,k+1)-vym(i,j,k-1))/dzm(k)/2
        vzdzm(i,j,k)=(vzm(i,j,k+1)-vzm(i,j,k-1))/dzm(k)/2
      enddo
      enddo
      enddo


      return
      end
