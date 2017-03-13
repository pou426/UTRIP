c======================================================================|
      subroutine psolv(work,i1,i2,dmat,cmat,margin,ix,jx,kx)
c======================================================================|
c
c NAME  psolv
c
c PURPOSE
c     solve x=M^-1 y, here M=LDU is precoditioner
c
c OUTPUTS
c     work(ix,jx,kx,7): [double] dimensions for vectors
c
c INPUTS
c     ix,jx,kx: [integer] dimension size
c     margin: [integer] size of boundary margins
c     cmat(ix,jx,kx,7): [double] sparse coeffcient matrix, "A"
c     dmat(ix,jx,kx): [double] ILU decomposed diagonal matrix
c
c HISTORY
c    written 2004-3-26 K. Nakamura
c    vectorized 2006-2-28 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension work(ix,jx,kx,7),dmat(ix,jx,kx)
      dimension cmat(ix,jx,kx,7)
      dimension zz(ix,jx,kx)
      integer   i1,i2
c----------------------------------------------------------------------|
      lh = ix-margin + jx-margin + kx-margin
      ll =  1+margin +  1+margin +  1+margin
c----------------------------------------------------------------------|
c     backward substitution
c     z=L^-1 y
      
      do k=1,kx
      do j=1,jx
      do i=1,ix
         zz(i,j,k)=0.
      enddo
      enddo
      enddo

c      do k=1+margin,kx-margin
c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         zz(i,j,k)=work(i,j,k,i2)
c     &        -cmat(i,j,k,2)*zz(i-1,j,k)-cmat(i,j,k,4)*zz(i,j-1,k)
c     &        -cmat(i,j,k,6)*zz(i,j,k-1)
c         zz(i,j,k)=zz(i,j,k)*dmat(i,j,k)
c      enddo
c      enddo
c      enddo
      do l=ll,lh
       ml = max( ll-( 1+margin), l-(ix-margin))
       mh = min( lh-(ix-margin), l-( 1+margin))
      do m=ml,mh
       kl = max( 1+margin, m-(jx-margin))
       kh = min(kx-margin, m-( 1+margin))
      do k=kl,kh
         i=l-m
         j=m-k
         zz(i,j,k)=dmat(i,j,k)*(
     &        work(i,j,k,i2)
     &        -cmat(i,j,k,2)*zz(i-1,  j,  k)
     &        -cmat(i,j,k,4)*zz(  i,j-1,  k)
     &        -cmat(i,j,k,6)*zz(  i,  j,k-1)
     &        )
      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|
c     forward substitution
c     x=(DU)^-1 z

      do k=1,kx   
      do j=1,jx
      do i=1,ix
         work(i,j,k,i1)=0.
      enddo
      enddo
      enddo

c      
c      do k=kx-margin,1+margin,-1
c      do j=jx-margin,1+margin,-1
c      do i=ix-margin,1+margin,-1
c         work(i,j,k,i1)=
c     &        cmat(i,j,k,7)*work(i,j,k+1,i1) + 
c     &        cmat(i,j,k,5)*work(i,j+1,k,i1) +
c     &        cmat(i,j,k,3)*work(i+1,j,k,i1)
c         work(i,j,k,i1)=zz(i,j,k)-dmat(i,j,k)*work(i,j,k,i1)
c      enddo
c      enddo
c      enddo
      do l=lh,ll,-1
       ml=max(ll-( 1+margin),l-(ix-margin))
       mh=min(lh-(ix-margin),l-( 1+margin))
      do m=mh,ml,-1
       kl=max( 1+margin, m-(jx-margin))
       kh=min(kx-margin, m-( 1+margin))
      do k=kh,kl,-1
         i=l-m
         j=m-k
         work(i,j,k,i1)=zz(i,j,k)-dmat(i,j,k)*(
     &         cmat(i,j,k,7)*work(i,j,k+1,i1)
     &        +cmat(i,j,k,5)*work(i,j+1,k,i1)
     &        +cmat(i,j,k,3)*work(i+1,j,k,i1) )
      enddo
      enddo
      enddo
      
c----------------------------------------------------------------------|

      return
      end
