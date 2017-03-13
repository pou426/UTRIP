c======================================================================|
      subroutine gtnorm(sum,res,margin,ix,jx)
c======================================================================|
c
c NAME  gtnorm
c
c PURPOSE
c    calculate norm of residue vector
c
c OUTPUTS
c    sum: [double] square of norm of residue vector
c    
c INPUTS
c    ix,jx: [integer] dimension size
c    res(ix,jx) : [double] residue
c    margin: [integer] size of boundary margins
c
c HISTORY
c    written 2002-3-1 T. Yokoyama 
c 
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension res(ix,jx)
c----------------------------------------------------------------------|
       sum=0.0e0
       do j=1+margin,jx-margin
       do i=1+margin,ix-margin
          sum=sum+res(i,j)**2
       enddo
       enddo


       return
       end
