c======================================================================|
      subroutine iludcmp(cmat,dmat,margin,ix,jx,kx)
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
c     ix,jx,kx: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,jx,kx,7): [double] sparse coefficient matrix, "A"
c 
c HISTORY
c    written 2004-3-26 K. Nakamura
c    vectorized 2006-2-28 K. Nakamura
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7)
      dimension dmat(ix,jx,kx)
c----------------------------------------------------------------------|
      lh = ix-margin + jx-margin + kx-margin
      ll =  1+margin +  1+margin +  1+margin
c----------------------------------------------------------------------|
      do k=1,kx
      do j=1,jx
      do i=1,ix
         dmat(i,j,k)=0.
      enddo
      enddo
      enddo

c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         dmat(i,j,k)=1./(cmat(i,j,k,1)
c     &        -cmat(i,j,k,2)*dmat(i-1,j,k)*cmat(i-1,j,k,3)
c     &        -cmat(i,j,k,4)*dmat(i,j-1,k)*cmat(i,j-1,k,5)
c     &        -cmat(i,j,k,6)*dmat(i,j,k-1)*cmat(i,j,k-1,7))
c      enddo
c      enddo
c      enddo
      do l=ll,lh
       ml=max(ll-( 1+margin),l-(ix-margin))
       mh=min(lh-(ix-margin),l-( 1+margin))
      do m=ml,mh
       kl=max( 1+margin, m-(jx-margin))
       kh=min(kx-margin, m-( 1+margin))
      do k=kl,kh
         i=l-m
         j=m-k
         dmat(i,j,k)=1./( cmat(i,j,k,1)
     &        -cmat(i,j,k,2)*dmat(i-1,j,k)*cmat(i-1,j,k,3)
     &        -cmat(i,j,k,4)*dmat(i,j-1,k)*cmat(i,j-1,k,5)
     &        -cmat(i,j,k,6)*dmat(i,j,k-1)*cmat(i,j,k-1,7) )
      enddo
      enddo
      enddo

      return
      end
