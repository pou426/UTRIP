c======================================================================|
      subroutine bdfdvy(mbnd,margin,bx,by,y,dx,ix,jx)
c======================================================================|
c
c NAME  bdfdvy
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
      dimension y(jx),dx(ix)
c----------------------------------------------------------------------|

      if (mbnd.eq.0) then 
        jbnd=1+margin
        do j=1,margin 
        do i=2,ix-1
          dbxdx=(bx(i+1,jbnd)-bx(i-1,jbnd))/dx(i)/2
          by(i,jbnd-j) = by(i,jbnd)-dbxdx*(y(jbnd-j)-y(jbnd))
        enddo
          i=1
          dbxdx=(bx(i+1,jbnd)-bx(i,jbnd))/(dx(i+1)+dx(i))*2
          by(i,jbnd-j) = by(jbnd,j)-dbxdx*(y(jbnd-j)-y(jbnd))
          i=ix
          dbxdx=(bx(i,jbnd)-bx(i-1,jbnd))/(dx(i)+dx(i-1))*2
          by(i,jbnd-j) = by(i,jbnd)-dbxdx*(y(jbnd-j)-y(jbnd))
        enddo
      else
        jbnd=jx-margin
        do j=1,margin
        do i=2,ix-1
          dbxdx=(bx(i+1,jbnd)-bx(i-1,jbnd))/dx(i)/2
          by(i,jbnd+j) = by(i,jbnd)-dbxdx*(y(jbnd+j)-y(jbnd))
        enddo
          i=1
          dbxdx=(bx(i+1,jbnd)-bx(i,jbnd))/(dx(i+1)+dx(i))*2
          by(i,jbnd+j) = by(i,jbnd)-dbxdx*(y(jbnd+j)-y(jbnd))
          i=ix
          dbxdx=(bx(i,jbnd)-bx(i-1,jbnd))/(dx(i)+dx(i-1))*2
          by(i,jbnd+j) = by(i,jbnd)-dbxdx*(y(jbnd+j)-y(jbnd))
        enddo
      endif
      
      return
      end
