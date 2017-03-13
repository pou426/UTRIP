c======================================================================|
      subroutine bnd(pot,x,y,margin,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pot(ix,jx)
      dimension x(ix),y(jx)
c----------------------------------------------------------------------|      

      potx0=0.d0
      potx1=1.d0
      dpotdx=(potx1-potx0)/(x(ix-margin+1)-x(margin))

c-- y0-boundary
      do i=1,ix
        do j=1,margin
          pot(i,j)=potx0+(x(i)-x(margin))*dpotdx*2.d0
        enddo
      enddo

c-- y1-boundary
      do i=1,ix
        do j=jx-margin+1,jx
          pot(i,j)=potx0+(x(i)-x(margin))*dpotdx
        enddo
      enddo

c-- x1-boundary
      do j=1,jx
        do i=ix-margin+1,ix
          pot(i,j)=potx1*2.d0-(y(j)-y(margin))*dpotdx*2.d0
          if (pot(i,j).lt.potx1) pot(i,j)=potx1
        enddo
      enddo

c-- x0-boundary
      do j=1,jx
        do i=1,margin
          pot(i,j)=potx0
        enddo
      enddo

      return
      end
