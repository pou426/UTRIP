c======================================================================|
      subroutine cipdxsrc(dadx,dady,da,dah,u,v,dt,dx,dy,ix,jx)
c======================================================================|
c
c NAME  cipdxsrc
c
c PURPOSE
c    advance no-advective (source-term) phase of CIP method for
c    physical variable gradients
c
c INPUTS & OUTPUTS
c    dadx(ix,jx): [double] physical variable gradient
c    dady(ix,jx): [double] physical variable gradient
c
c OUTPUTS
c    None
c
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c
c    da(ix,jx): [double] physical variable
c    dah(ix,jx): [double] physical variable after non-advective phase
c    u(ix,jx): [double] advection velocity
c    v(ix,jx): [double] advection velocity
c    dt: [double] delta time
c    dx(ix),dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on T. Kudoh's code
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension dadx(ix,jx),dady(ix,jx),u(ix,jx),v(ix,jx)
      dimension da(ix,jx),dah(ix,jx)
      dimension dadxn(ix,jx),dadyn(ix,jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
        dadxn(i,j)=0.d0
        dadyn(i,j)=0.d0
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
        dadxn(i,j)= dadx(i,j)
     &         +((dah(i+1,j)-da(i+1,j))-(dah(i-1,j)-da(i-1,j)))/2/dx(i)
     &            - dadx(i,j)*(u(i+1,j)-u(i-1,j))*dt/dx(i)/2 
     &            - dady(i,j)*(v(i+1,j)-v(i-1,j))*dt/dx(i)/2 
        dadyn(i,j)= dady(i,j)
     &         +((dah(i,j+1)-da(i,j+1))-(dah(i,j-1)-da(i,j-1)))/2/dy(j)
     &            - dadx(i,j)*(u(i,j+1)-u(i,j-1))*dt/dy(j)/2 
     &            - dady(i,j)*(v(i,j+1)-v(i,j-1))*dt/dy(j)/2 
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
        dadx(i,j)=dadxn(i,j)
        dady(i,j)=dadyn(i,j)
      enddo
      enddo

      return
      end
