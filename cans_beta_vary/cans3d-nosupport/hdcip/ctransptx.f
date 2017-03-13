c======================================================================|
      subroutine ctransptx(bxmh,vxcy,vxcz,bxcy,bxcz,vycx,vzcx,bycx,bzcx
     &                         ,dt,dy,dz,ix,jx,kx)
c======================================================================|
c
c NAME  ctransptx
c
c PURPOSE
c    Constrained Transport (CT) method for induction equation
c    satisfying div B = 0.
c
c INPUTS & OUTPUTS
c    bxh(ix,jx,kx): [double] magnetic field
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    vym(ix,jx,kx): [double] velocity
c    bym(ix,jx,kx): [double] magnetic field
c    vxm(ix,jx,kx) : [double] velocity along the x-cordinate
c    bxm(ix,jx,kx) : [double] magnetic field
c    dx(ix) : [double] grid spacing
c    dt: [double] delta time
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dy(jx),dz(kx)

      dimension bxmh(ix,jx,kx)
      dimension eyc(ix,jx,kx),ezc(ix,jx,kx)
      dimension vxcy(ix,jx,kx),vxcz(ix,jx,kx)
      dimension bxcy(ix,jx,kx),bxcz(ix,jx,kx)
      dimension vycx(ix,jx,kx),vzcx(ix,jx,kx)
      dimension bycx(ix,jx,kx),bzcx(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        eyc(i,j,k)=-(vzcx(i,j,k)*bxcz(i,j,k)-vxcz(i,j,k)*bzcx(i,j,k))
        ezc(i,j,k)=-(vxcy(i,j,k)*bycx(i,j,k)-vycx(i,j,k)*bxcy(i,j,k))
      enddo
      enddo
      enddo

      do k=2,kx-1
      do j=2,jx-1
      do i=1,ix-1
        bxmh(i,j,k)=bxmh(i,j,k)-dt/dy(j)*(ezc(i,j,k)-ezc(i,j-1,k))
     &                         +dt/dz(k)*(eyc(i,j,k)-eyc(i,j,k-1))
      enddo
      enddo
      enddo

      return
      end
