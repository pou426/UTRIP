c======================================================================|
      subroutine sorbr(te,res,mrb,omsor,cmat,src,margin,ix,jx,kx)
c======================================================================|
c
c NAME  sorbr
c 
c PURPOSE
c    solve matrix inversion using SOR red & black method
c    
c INPUTS & OUTPUTS
c    te(ix,jx,kx): [double] temperature 
c        
c OUTPUTS
c    res(ix,jx,kx) : [double] residue
c    
c INPUTS
c    ix,jx,kx: [integer] dimension size
c    mrb: [integer] mrb=1 for red step, mrb=2 for black step
c    omsor: [double] over-relaxation parameter (1<omsor<2)
c    cmat(ix,jx,kx,7): [double] coefficient matrix of heat conduction eq
c    src(ix,jx,kx): [double] source vector of heat conduction eq
c    margin: [integer] size of boundary margins
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
c      
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),src(ix,jx,kx)
      dimension res(ix,jx,kx)
c----------------------------------------------------------------------|

c    if mrb=1, (i,j,k)=(e,o,o)
c    if mrb=2, (i,j,k)=(o,o,o)

       do k=margin+1,kx-margin,2
       do j=margin+1,jx-margin,2
       do i=margin+(3-mrb),ix-margin,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo

c    if mrb=1, (i,j,k)=(o,o,e)
c    if mrb=2, (i,j,k)=(e,o,e)

       do k=margin+2,kx-margin,2
       do j=margin+1,jx-margin,2
       do i=margin+mrb,ix-margin,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo


c    if mrb=1, (i,j,k)=(o,e,o)
c    if mrb=2, (i,j,k)=(e,e,o)

       do k=margin+1,kx-margin,2
       do j=margin+2,jx-margin,2
       do i=margin+mrb,ix-margin,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo

c    if mrb=1, (i,j,k)=(e,e,e)
c    if mrb=2, (i,j,k)=(o,e,e)

       do k=margin+2,kx-margin,2
       do j=margin+2,jx-margin,2
       do i=margin+(3-mrb),ix-margin,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo


      return
      end
