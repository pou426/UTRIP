c======================================================================|
      subroutine gtnorm2(sum,res,margar,ix,jx,kx)
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
c    margar: [integer] size of boundary margar
c
c HISTORY
c    written 2002-3-1 T. Yokoyama 
c 
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension margar(3)
      dimension res(ix,jx,kx)
c----------------------------------------------------------------------|
       sum=0.0d0
       do k=1+margar(3),kx-margar(3)
       do j=1+margar(2),jx-margar(2)
       do i=1+margar(1),ix-margar(1)
          sum=sum+res(i,j,k)**2
       enddo
       enddo
       enddo

       return
       end
