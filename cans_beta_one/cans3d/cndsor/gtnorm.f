c======================================================================|
      subroutine gtnorm(sum,res,margin,ix,jx,kx)
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
c    ix,jx,kx: [integer] dimension size
c    res(ix,jx,kx) : [double] residue
c    margin: [integer] size of boundary margins
c
c HISTORY
c    written 2002-3-1 T. Yokoyama 
c 
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension res(ix,jx,kx)
c----------------------------------------------------------------------|
       sum=0.0e0
       do k=1+margin,kx-margin
       do j=1+margin,jx-margin
       do i=1+margin,ix-margin
          sum=sum+res(i,j,k)**2
       enddo
       enddo
       enddo



       return
       end
