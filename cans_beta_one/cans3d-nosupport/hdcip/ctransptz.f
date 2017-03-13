c======================================================================|
      subroutine ctransptz(bzmh,vzcx,vzcy,bzcx,bzcy,vxcz,bxcz,vycz,bycz
     &                         ,dt,dx,dy,ix,jx,kx)
c======================================================================|
c
c NAME  ctransptz
c
c PURPOSE
c    Constrained Transport (CT) method for induction equation
c    satisfying div B = 0.
c
c INPUTS & OUTPUTS
c    bzmh(ix,jx,kx): [double] magnetic field
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
c    NOTE: T. Yokoyama's code is based on T. Kudoh's code
c    written 2003-6-7 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|

      implicit real*8 (a-h,o-z)

      dimension dx(ix),dy(jx)

      dimension bzmh(ix,jx,kx)
      dimension exc(ix,jx,kx),eyc(ix,jx,kx)
      dimension vzcx(ix,jx,kx),vzcy(ix,jx,kx)
      dimension bzcx(ix,jx,kx),bzcy(ix,jx,kx)
      dimension vxcz(ix,jx,kx),vycz(ix,jx,kx)
      dimension bxcz(ix,jx,kx),bycz(ix,jx,kx) 
c----------------------------------------------------------------------|

      do k=1,kx-1
      do j=1,jx-1
      do i=1,ix-1
        exc(i,j,k)=-(vycz(i,j,k)*bzcy(i,j,k)-vzcy(i,j,k)*bycz(i,j,k))
        eyc(i,j,k)=-(vzcx(i,j,k)*bxcz(i,j,k)-vxcz(i,j,k)*bzcx(i,j,k))
      enddo
      enddo
      enddo

      do k=1,kx-1
      do j=2,jx-1
      do i=2,ix-1
        bzmh(i,j,k)=bzmh(i,j,k)-dt/dx(i)*(eyc(i,j,k)-eyc(i-1,j,k))
     &                         +dt/dy(j)*(exc(i,j,k)-exc(i,j-1,k))
      enddo
      enddo
      enddo

      return
      end
