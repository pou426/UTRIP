c======================================================================|
      subroutine residue(anorm,res,margin,ix)
c======================================================================|
c
c NAME  residue
c
c PURPOSE
c    calculate norm of residue vector
c
c OUTPUTS
c    anorm: [double] norm of residue vector
c
c INPUTS
c    ix: [integer] dimension size
c    res(ix) : [double] residue
c    margin: [integer] size of boundary margins
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension res(ix)
c----------------------------------------------------------------------|
       sum=0.0e0
       do i=1+margin,ix-margin
          sum=sum+res(i)**2
       enddo
       anorm=sqrt(sum)



       return
       end
