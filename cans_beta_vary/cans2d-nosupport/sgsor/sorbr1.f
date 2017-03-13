c======================================================================|
      subroutine sorbr1(te,res,mrb,omsor,cmat,src,margin,ix,jx)
c======================================================================|
c
c NAME  sorbr
c
c PURPOSE
c    solve matrix inversion using SOR red & black method
c
c INPUTS & OUTPUTS
c    te(ix,jx): [double] temperature
c    
c OUTPUTS
c    res(ix,jx) : [double] residue
c    
c INPUTS
c    ix,jx: [integer] dimension size
c    mrb: [integer] mrb=1 for red step, mrb=2 for black step
c    omsor: [double] over-relaxation parameter (1<omsor<2)
c    cmat(ix,jx,5): [double] coefficient matrix of heat conduction eq
c    src(ix,jx): [double] source vector of heat conduction eq
c    margin: [integer] size of boundary margins
c     
c HISTORY
c    written 2002-3-1 T. Yokoyama
c      
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension te(ix,jx)
      dimension cmat(ix,jx,5),src(ix,jx)
      dimension res(ix,jx)
c----------------------------------------------------------------------|

c    if mrb=1, (i,j)=(e,o)
c    if mrb=2, (i,j)=(o,o)

       do j=margin+1,jx-margin,2
       do i=margin+(3-mrb),ix-margin,2
          res(i,j)=  cmat(i,j,1)*te(i,j)-src(i,j)
     &              +cmat(i,j,2)*te(i-1,j)+cmat(i,j,3)*te(i+1,j)
     &              +cmat(i,j,4)*te(i,j-1)+cmat(i,j,5)*te(i,j+1)
          tmp=te(i,j)
          te(i,j)=te(i,j)-omsor*res(i,j)/cmat(i,j,1)
c         if (i.eq.30.and.j.eq.30) 
c    &      write(6,*) i,j,res(i,j),cmat(i,j,1)*te(i,j),src(i,j)
c    &              ,cmat(i,j,2)*te(i-1,j),cmat(i,j,3)*te(i+1,j)
c    &              ,cmat(i,j,4)*te(i,j-1),cmat(i,j,5)*te(i,j+1)
       enddo
       enddo

c    if mrb=1, (i,j)=(o,e)
c    if mrb=2, (i,j)=(e,e)

       do j=margin+2,jx-margin,2
       do i=margin+mrb,ix-margin,2
          res(i,j)=  cmat(i,j,1)*te(i,j)-src(i,j)
     &              +cmat(i,j,2)*te(i-1,j)+cmat(i,j,3)*te(i+1,j)
     &              +cmat(i,j,4)*te(i,j-1)+cmat(i,j,5)*te(i,j+1)
          te(i,j)=te(i,j)-omsor*res(i,j)/cmat(i,j,1)
       enddo
       enddo

      return
      end
