c======================================================================|
      subroutine sorbr(te,res,mrb,omsor,cmat,src,margin,ix)
c======================================================================|
c
c NAME  sorbr
c
c PURPOSE
c    solve matrix inversion using SOR red & black method
c
c INPUTS & OUTPUTS
c    te(ix): [double] temperature
c
c OUTPUTS
c    res(ix) : [double] residue
c
c INPUTS
c    ix: [integer] dimension size
c    mrb: [integer] mrb=1 for red step, mrb=2 for black step
c    omsor: [double] over-relaxation parameter (1<omsor<2)
c    cmat(ix,3): [double] coefficient matrix of heat conduction eq
c    src(ix): [double] source vector of heat conduction eq
c    margin: [integer] size of boundary margins
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension te(ix)
      dimension cmat(ix,3),src(ix)
      dimension res(ix)
c----------------------------------------------------------------------|

c    if mrb=1, (i)=(o)
c    if mrb=2, (i)=(e)

       do i=margin+mrb,ix-margin,2
          res(i)=  cmat(i,1)*te(i)-src(i)
     &            +cmat(i,2)*te(i-1)+cmat(i,3)*te(i+1)
          te(i)=te(i)-omsor*res(i)/cmat(i,1)
       enddo


      return
      end
