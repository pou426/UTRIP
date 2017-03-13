c======================================================================|
      subroutine cipdxsrc(dadx,dady,dadz,da,dah,u,v,w,dt
     &                        ,dx,dy,dz,ix,jx,kx)
c======================================================================|
c
c NAME  cipdxsrc
c
c PURPOSE
c    advance no-advective (source-term) phase of CIP method for
c    physical variable gradients
c
c INPUTS & OUTPUTS
c    dadx(ix,jx,kx): [double] physical variable gradient
c    dady(ix,jx,kx): [double] physical variable gradient
c    dadz(ix,jx,kx): [double] physical variable gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    da(ix,jx,kx): [double] physical variable
c    dah(ix,jx,kx): [double] physical variable after non-advective phase
c    u(ix,jx,kx): [double] advection velocity
c    v(ix,jx,kx): [double] advection velocity
c    w(ix,jx,kx): [double] advection velocity
c    dt: [double] delta time
c    dx(ix),dy(jx),dz(kx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2003-6-1 K. Takahashi based on T. Yokoyama's code
c
c----------------------------------------------------------------------|
      implicit real*8 (a-h,o-z)
      dimension dx(ix),dy(jx),dz(kx)
      dimension dadx(ix,jx,kx),dady(ix,jx,kx),dadz(ix,jx,kx)
      dimension u(ix,jx,kx),v(ix,jx,kx),w(ix,jx,kx)
      dimension da(ix,jx,kx),dah(ix,jx,kx)
      dimension dadxn(ix,jx,kx),dadyn(ix,jx,kx),dadzn(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
        dadxn(i,j,k)=0.d0
        dadyn(i,j,k)=0.d0
        dadzn(i,j,k)=0.d0
      enddo
      enddo
      enddo

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
        dadxn(i,j,k)= dadx(i,j,k)
     &              +((dah(i+1,j,k)-da(i+1,j,k))
     &               -(dah(i-1,j,k)-da(i-1,j,k)))/2/dx(i)
     &              - dadx(i,j,k)*(u(i+1,j,k)-u(i-1,j,k))*dt/dx(i)/2 
     &              - dady(i,j,k)*(v(i+1,j,k)-v(i-1,j,k))*dt/dx(i)/2 
     &              - dadz(i,j,k)*(w(i+1,j,k)-w(i-1,j,k))*dt/dx(i)/2 
        dadyn(i,j,k)= dady(i,j,k)
     &              +((dah(i,j+1,k)-da(i,j+1,k))
     &               -(dah(i,j-1,k)-da(i,j-1,k)))/2/dy(j)
     &              - dadx(i,j,k)*(u(i,j+1,k)-u(i,j-1,k))*dt/dy(j)/2 
     &              - dady(i,j,k)*(v(i,j+1,k)-v(i,j-1,k))*dt/dy(j)/2 
     &              - dadz(i,j,k)*(w(i,j+1,k)-w(i,j-1,k))*dt/dy(j)/2
        dadzn(i,j,k)= dadz(i,j,k)
     &              +((dah(i,j,k+1)-da(i,j,k+1))
     &               -(dah(i,j,k-1)-da(i,j,k-1)))/2/dz(k)
     &              - dadx(i,j,k)*(u(i,j,k+1)-u(i,j,k-1))*dt/dz(k)/2 
     &              - dady(i,j,k)*(v(i,j,k+1)-v(i,j,k-1))*dt/dz(k)/2 
     &              - dadz(i,j,k)*(w(i,j,k+1)-w(i,j,k-1))*dt/dz(k)/2
      enddo
      enddo
      enddo

      do k=2,kx-1
      do j=2,jx-1
      do i=2,ix-1
        dadx(i,j,k)=dadxn(i,j,k)
        dady(i,j,k)=dadyn(i,j,k)
        dadz(i,j,k)=dadzn(i,j,k)
      enddo
      enddo
      enddo

      return
      end
