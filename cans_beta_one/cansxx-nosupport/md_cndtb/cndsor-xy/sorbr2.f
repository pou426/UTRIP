c======================================================================|
      subroutine sorbr2(te,res,mrb,omsor,cmat,src,ix,jx,kx,mfdim,margar)
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
c    mrb: [integer] mrb=0 for red step, mrb=1 for black step
c    omsor: [double] over-relaxation parameter (1<omsor<2)
c    cmat(ix,jx,kx,7): [double] coefficient matrix of heat conduction eq
c    src(ix,jx,kx): [double] source vector of heat conduction eq
c    margar: [integer] size of boundary margar
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
c      
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension mfdim(3),margar(3)
      dimension te(ix,jx,kx)
      dimension cmat(ix,jx,kx,7),src(ix,jx,kx)
      dimension res(ix,jx,kx)
c----------------------------------------------------------------------|
      i0=margar(1)+1
      i1=ix-margar(1)
      j0=margar(2)+1
      j1=jx-margar(2)
      k0=margar(3)+1
      k1=kx-margar(3)

c    if mrb=0, (i,j,k)=(e,o,o)
c    if mrb=1, (i,j,k)=(o,o,o)

       do k=k0,k1,2
       do j=j0,j1,2
       do i=i0+(1-mrb),i1,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo

c    if mrb=0, (i,j,k)=(o,o,e)
c    if mrb=1, (i,j,k)=(e,o,e)

       do k=k0+1,k1,2
       do j=j0,j1,2
       do i=i0+mrb,i1,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo


c    if mrb=0, (i,j,k)=(o,e,o)
c    if mrb=1, (i,j,k)=(e,e,o)

       do k=k0,k1,2
       do j=j0+1,j1,2
       do i=i0+mrb,i1,2
         res(i,j,k)= cmat(i,j,k,1)*te(i,j,k)-src(i,j,k)
     &              +cmat(i,j,k,2)*te(i-1,j,k)+cmat(i,j,k,3)*te(i+1,j,k)
     &              +cmat(i,j,k,4)*te(i,j-1,k)+cmat(i,j,k,5)*te(i,j+1,k)
     &              +cmat(i,j,k,6)*te(i,j,k-1)+cmat(i,j,k,7)*te(i,j,k+1)
          te(i,j,k)=te(i,j,k)-omsor*res(i,j,k)/cmat(i,j,k,1)
       enddo
       enddo
       enddo

c    if mrb=0, (i,j,k)=(e,e,e)
c    if mrb=1, (i,j,k)=(o,e,e)

       do k=k0+1,k1,2
       do j=j0+1,j1,2
       do i=i0+(1-mrb),i1,2
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
