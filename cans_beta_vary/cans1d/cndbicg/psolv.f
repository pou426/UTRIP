c======================================================================|
      subroutine psolv(work,i1,i2,dmat,cmat,margin,ix)
c======================================================================|
c
c NAME  psolv
c
c PURPOSE
c     solve x=M^-1 y, here M=LDU is preconditioner
c
c OUTPUTS
c     work(i,i1): [double] vector x
c
c INPUTS
c     ix: [integer] dimension size
c     margin: [integer] size of boundary margin
c     work(i,i2): [double] vector y
c     i1: [integer] index for vector x
c     i2: [integer] index for vector y
c     dmat: [double] ILU decomposed diagonal matrix
c     cmat: [double] coefficient matrix
c
c HISTORY
c    written 2002-3-23 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension work(ix,1:7),dmat(ix)
      dimension cmat(ix,1:3)
      dimension zz(ix)
      integer i1,i2
C----------------------------------------------------------------------|
c     LDUx=y => Lz=y & DUx=z

c     backward substitution
c     z=L^-1 y
      do i=1,margin
        zz(i)=0.
      enddo
      do i=1+margin,ix-margin
         zz(i) = work(i,i2) - cmat(i,2)*zz(i-1)
         zz(i) = zz(i)*dmat(i)
      enddo


c     forward substitution
c     x=(DU)^-1 z
      do i=ix+1-margin,ix
        work(i,i1)=0.
      enddo
      do i=ix-margin,1+margin,-1
         work(i,i1) = cmat(i,3)*work(i+1,i1)
         work(i,i1) = zz(i) - dmat(i)*work(i,i1)
      enddo


      return
      end
