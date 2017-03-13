c======================================================================|
      subroutine prtote(te,ro,pr,gm,ix,jx)
c======================================================================|
c
c NAME  prtote
c
c PURPOSE
c    calculation of temperature
c
c OUTPUTS
c    te(ix,jx): [double] temperature
c
c INPUTS
c    ro(ix,jx): [double] density
c    pr(ix,jx): [double] pressure
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
            te(i,j)=gm*pr(i,j)/ro(i,j)
         enddo
         enddo

         return
         end
