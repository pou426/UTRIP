c======================================================================|
      subroutine pertub(vy,x,ix,y,jx,margin)
c======================================================================|
      implicit real*8 (a-h,o-z)
      dimension x(ix), y(jx)
      dimension vy(ix,jx)

c----------------------------------------------------------------------|
c     perturbation
c----------------------------------------------------------------------|
      pi = acos(-1.0d0)
      marginx=5

      amp=1.e-7
      xmin=x(marginx+1)
      xmax=x(ix-marginx)
      ymin=y(margin+1)
      ymax=y(jx-margin)
      a1=xmax-xmin
      a2=ymax-ymin

      do jp=1,1
      do ip=1,1
      do j=1,jx
      do i=1,ix
       xx=x(i)-xmin
       yy=y(j)-ymin
       vy(i,j) =  vy(i,j)
     &    +amp*sin(2*ip*pi*xx/a1)*sin(jp*pi*yy/a2)

      enddo
      enddo
      enddo
      enddo

      return
      end



