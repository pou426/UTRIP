c======================================================================|
      subroutine bbtocx_s(cx,bz,dy,x,y,ix,jx)
c======================================================================|
c
c NAME  bbtocx_s
c
c PURPOSE
c    calculate current density
c        * r-component in cylindrical coordinate
c
c OUTPUTS
c    cx(ix,jx): [double] current density
c
c INPUTS
c    bz(ix,jx): [double] magnetic field
c    x(ix),y(jx) : [double] coordinate
c    dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension x(ix),y(jx)
      dimension dy(jx)
      dimension bz(ix,jx)
      dimension cx(ix,jx)
c----------------------------------------------------------------------|

         do j=2,jx-1
         do i=2,ix-1
           cx(i,j) =  (x(i)*sin(y(j+1))*bz(i,j+1)
     &                -x(i)*sin(y(j-1))*bz(i,j-1))/dy(j)/2
     &                /x(i)**2/sin(y(j))
         enddo
         enddo

          return
          end
