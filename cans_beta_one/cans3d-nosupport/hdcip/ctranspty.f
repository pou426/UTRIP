c======================================================================|
      subroutine ctranspty(bymh,vycx,vycz,bycx,bycz,vxcy,vzcy,bxcy,bzcy
     &                         ,dt,dx,dz,ix,jx,kx)
c======================================================================|
c
c NAME  ctranspty
c
c PURPOSE
c    Constrained Transport (CT) method for induction equation
c    satisfying div B = 0.
c
c INPUTS & OUTPUTS
c    byh(ix,jx,kx): [double] magnetic field
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

      dimension dx(ix),dz(kx)

      dimension bymh(ix,jx,kx)
      dimension exc(ix,jx,kx),ezc(ix,jx,kx)
      dimension vycx(ix,jx,kx),vycz(ix,jx,kx)
      dimension bycx(ix,jx,kx),bycz(ix,jx,kx)
      dimension vxcy(ix,jx,kx),vzcy(ix,jx,kx)
      dimension bxcy(ix,jx,kx),bzcy(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        exc(i,j,k)=-(vycz(i,j,k)*bzcy(i,j,k)-vzcy(i,j,k)*bycz(i,j,k))
        ezc(i,j,k)=-(vxcy(i,j,k)*bycx(i,j,k)-vycx(i,j,k)*bxcy(i,j,k))
      enddo
      enddo
      enddo

      do k=2,kx-1
      do j=1,jx-1
      do i=2,ix-1
        bymh(i,j,k)=bymh(i,j,k)-dt/dz(k)*(exc(i,j,k)-exc(i,j,k-1))
     &                         +dt/dx(i)*(ezc(i,j,k)-ezc(i-1,j,k))
      enddo
      enddo
      enddo

      return
      end
