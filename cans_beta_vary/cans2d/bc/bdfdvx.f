c======================================================================|
      subroutine bdfdvx(mbnd,margin,bx,by,x,dy,ix,jx)
c======================================================================|
c
c NAME  bdfdvx
c
c PURPOSE
c    apply free boundary condition
c    sutisfying divB=0 condition
c
c INPUTS & OUTPUTS
c    bx(ix,jx): [double] magnetic field
c    by(ix,jx): [double] magnetic field
c
c OUTPUTS
c    None
c
c INPUTS
c    ix,jx: [integer] dimension size
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    mbnd: [integer] If mbnd=0, smaller 'i' side. 
c                    If mbnd=1, larger  'i' side.
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    bug fixed 2006-7-8 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension bx(ix,jx),by(ix,jx)
      dimension x(ix),dy(jx)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        ibnd=1+margin
        do i=1,margin 
        do j=2,jx-1
          dbydy=(by(ibnd,j+1)-by(ibnd,j-1))/dy(j)/2
          bx(ibnd-i,j) = bx(ibnd,j)-dbydy*(x(ibnd-i)-x(ibnd))
        enddo
          j=1
          dbydy=(by(ibnd,j+1)-by(ibnd,j))/(dy(j+1)+dy(j))*2
          bx(ibnd-i,j) = bx(ibnd,j)-dbydy*(x(ibnd-i)-x(ibnd))
          j=jx
          dbydy=(by(ibnd,j)-by(ibnd,j-1))/(dy(j)+dy(j-1))*2
          bx(ibnd-i,j) = bx(ibnd,j)-dbydy*(x(ibnd-i)-x(ibnd))
        enddo
      else
        ibnd=ix-margin
        do i=1,margin
        do j=2,jx-1
          dbydy=(by(ibnd,j+1)-by(ibnd,j-1))/dy(j)/2
          bx(ibnd+i,j) = bx(ibnd,j)-dbydy*(x(ibnd+i)-x(ibnd))
        enddo
          j=1
          dbydy=(by(ibnd,j+1)-by(ibnd,j))/(dy(j+1)+dy(j))*2
          bx(ibnd+i,j) = bx(ibnd,j)-dbydy*(x(ibnd+i)-x(ibnd))
          j=jx
          dbydy=(by(ibnd,j)-by(ibnd,j-1))/(dy(j)+dy(j-1))*2
          bx(ibnd+i,j) = bx(ibnd,j)-dbydy*(x(ibnd+i)-x(ibnd))
        enddo
      endif
      
      return
      end
