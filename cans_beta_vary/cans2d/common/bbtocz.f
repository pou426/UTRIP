c======================================================================|
      subroutine bbtocz(cz,bx,by,dx,dy,ix,jx)
c======================================================================|
c
c NAME  bbtocz
c
c PURPOSE
c    calculate current density
c        * z-component
c
c OUTPUTS
c    cz(ix,jx): [double] current density
c
c INPUTS
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c    dx(ix),dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension cz(ix,jx)
c----------------------------------------------------------------------|

         do j=2,jx-1
         do i=2,ix-1
c          cz(i,j) =  (by(i+1,j)-by(i-1,j))/dx(i)/2.
c    &               -(bx(i,j+1)-bx(i,j-1))/dy(j)/2.
           cz(i,j) =  (by(i+1,j)-by(i-1,j))/dx(i)/2.d0/2.d0
     &               +(by(i+1,j+1)-by(i-1,j+1))/dx(i)/2.d0/4.d0
     &               +(by(i+1,j-1)-by(i-1,j-1))/dx(i)/2.d0/4.d0
     &               -(bx(i,j+1)-bx(i,j-1))/dy(j)/2.d0/2.d0
     &               -(bx(i+1,j+1)-bx(i+1,j-1))/dy(j)/2.d0/4.d0
     &               -(bx(i-1,j+1)-bx(i-1,j-1))/dy(j)/2.d0/4.d0
         enddo
         enddo

          return
          end
