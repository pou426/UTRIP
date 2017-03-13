c======================================================================|
      subroutine iludcmp(cmat,dmat,margin,ix)
c======================================================================|
c
c NAME  iludcmp
c
c PURPOSE
c     decompose coefficient matrix to Imcomplete Upper-Lower matrix
c     Ax=y, A ~ M, M=LDU
c
c OUTPUTS
c     dmat(ix): [double] decomposed diagonal matrix, "D"
c
c INPUTS
c     ix: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,3): [double] sparse coefficient matrix, "A"
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,3)
      dimension dmat(ix)
C----------------------------------------------------------------------|
      do i=1,margin
        dmat(i)=0.
      enddo
      do i=1+margin,ix-margin
         dmat(i)=1./(cmat(i,1)-cmat(i,2)*dmat(i-1)*cmat(i-1,3))
      enddo

      return
      end
