c======================================================================|
      subroutine psolv(work,i1,i2,dmat,cmat,margin,ix,jx)
c======================================================================|
c
c NAME  psolv
c
c PURPOSE
c
c INPUTS & OUTPUTS
c
c OUTPUTS
c
c INPUTS
c
c HISTORY
c    written 2002-3-1 K. Nakamura
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension work(ix,jx,7),dmat(ix,jx)
      dimension cmat(ix,jx,5)
      dimension zz(ix,jx)
      integer   i1,i2
c----------------------------------------------------------------------|
      lh = ix-margin + jx-margin
      ll =  1+margin + 1+margin

c----------------------------------------------------------------------|
c     backward substitution
c     z=L^-1 y

      do j=1,jx
      do i=1,ix
         zz(i,j)=0.
      enddo
      enddo

c      do j=1+margin,jx-margin
c      do i=1+margin,ix-margin
c         tmp=work(i,j,i2)
c     &        -cmat(i,j,2)*zz(i-1,j)-cmat(i,j,4)*zz(i,j-1)
c         zz(i,j)=tmp*dmat(i,j)
c      enddo
c      enddo
      do l=ll,lh
         jl = max( 1+margin, l-(ix-margin))
         jh = min(jx-margin, l-(1+margin))
         do j=jl,jh
            i=l-j
            zz(i,j)=dmat(i,j)*( work(i,j,i2)
     &           -cmat(i,j,2)*zz(i-1,j)-cmat(i,j,4)*zz(i,j-1) )
         enddo
      enddo

c----------------------------------------------------------------------|
c     forward substitution
c     x=(DU)^-1 z

      do j=1,jx
      do i=1,ix
         work(i,j,i1)=0.
      enddo
      enddo

c      
c      do j=jx-margin,1+margin,-1
c      do i=ix-margin,1+margin,-1
c         tmp=cmat(i,j,5)*work(i,j+1,i1)+cmat(i,j,3)*work(i+1,j,i1)
c         work(i,j,i1)=zz(i,j)-dmat(i,j)*tmp
c      enddo
c      enddo
      do l=lh,ll,-1
         jl = max( 1+margin, l-(ix-margin))
         jh = min(jx-margin, l-(1+margin))
         do j=jh,jl,-1
            i=l-j
            work(i,j,i1)=zz(i,j)-dmat(i,j)*(
     &           cmat(i,j,5)*work(i,j+1,i1)+cmat(i,j,3)*work(i+1,j,i1) )
         enddo
      enddo

c----------------------------------------------------------------------|

      return
      end
