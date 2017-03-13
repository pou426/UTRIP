c======================================================================|
      subroutine tetopr(ro,pr,te,gm,ix,jx)
c======================================================================|
c
c NAME  tetopr
c
c PURPOSE
c    calculation of pressure from temperature
c
c OUTPUTS
c    pr(ix,jx): [double] pressure
c
c INPUTS
c    te(ix,jx): [double] temperature
c    ro(ix,jx): [double] density
c    gm: [double] polytropic index gamma
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)

      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension te(ix,jx)
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
        pr(i,j)=ro(i,j)*te(i,j)/gm
      enddo
      enddo

      return
      end
