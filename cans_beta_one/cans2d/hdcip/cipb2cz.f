c======================================================================|
      subroutine cipb2cz(czc,bxm,bym,dxm,dym,ix,jx)
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
c    bug-fixed 2006-3-12 T. Yokoyama (informed by H. Isobe)
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dym(jx)
      dimension bxm(ix,jx),bym(ix,jx)
      dimension czc(ix,jx)
c----------------------------------------------------------------------|

         do j=2,jx-1
         do i=2,ix-1
           czc(i,j) =  (bym(i+1,j)-bym(i,j))/dxm(i)
     &                -(bxm(i,j+1)-bxm(i,j))/dym(j)
         enddo
         enddo

          return
          end
