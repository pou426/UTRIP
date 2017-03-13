c======================================================================|
      subroutine iludcmp(cmat,dmat,margin,ix,jx)
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
c     ix,jx: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,jx,5): [double] sparse coefficient matrix, "A"
c 
c HISTORY
c    written 2002-3-23 K. Nakamura
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5)
      dimension dmat(ix,jx)
c----------------------------------------------------------------------|
      do j=1,jx
      do i=1,ix
         dmat(i,j)=0.
      enddo
      enddo

c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         dmat(i,j)=1./(cmat(i,j,1)
c     &        -cmat(i,j,2)*dmat(i-1,j)*cmat(i-1,j,3)
c     &        -cmat(i,j,4)*dmat(i,j-1)*cmat(i,j-1,5))
c      enddo
c      enddo
      ll =  1+margin +  1+margin
      lh = ix-margin + jx-margin
      do l=ll,lh
         jl = max( 1+margin ,l-(ix-margin))
         jh = min(jx-margin ,l- (1+margin))
         do j=jl,jh
            i=l-j
            dmat(i,j) = 1./( cmat(i,j,1)
     &           -cmat(i,j,2)*dmat(i-1,  j)*cmat(i-1,  j,3)
     &           -cmat(i,j,4)*dmat(i  ,j-1)*cmat(i  ,j-1,5) )
         enddo
      enddo

      return
      end
