c======================================================================|
      subroutine visten(visxx,visxy,visyy,vx,vy,dx,ix,dy,jx,visc)
c======================================================================|
c
c NAME  visten
c
c PURPOSE
c    calculate viscous tensor
c
c OUTPUTS
c     vis???(ix,jx): [double] viscous tensor
c
c INPUTS
c    v?(ix,jx): [double] velocity
c    dx(ix),dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2005/6/29 H.Isobe
c
c----------------------------------------------------------------------|
      implicit real*8 (a-h,o-z)
      dimension dx(ix),dy(jx)
      dimension vx(ix,jx),vy(ix,jx)
      dimension visxx(ix,jx),visxy(ix,jx),visyy(ix,jx)
      dimension visc(ix,jx)
c----------------------------------------------------------------------|
      twoth=2.d0/3.d0


         do j=3,jx-2
         do i=3,ix-2

           x1= dy(j)/2+dy(j+1)/2
           x2= dy(j)/2+dy(j+1)+dy(j+2)/2
           x3=-dy(j)/2-dy(j-1)/2
           x4=-dy(j)/2-dy(j-1)-dy(j-2)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=vx(i,j+1)-vx(i,j)
           y2=vx(i,j+2)-vx(i,j)
           y3=vx(i,j-1)-vx(i,j)
           y4=vx(i,j-2)-vx(i,j)

           w1=vy(i,j+1)-vy(i,j)
           w2=vy(i,j+2)-vy(i,j)
           w3=vy(i,j-1)-vy(i,j)
           w4=vy(i,j-2)-vy(i,j)

           a11=x1*(a3*(1-a3**3)*a2**3*(1-a2)
     &            -a2*(1-a2**3)*a3**3*(1-a3))
           a12=x1**2*(a3**2*(1-a3**2)*a2**3*(1-a2)
     &               -a2**2*(1-a2**2)*a3**3*(1-a3))
           a21=x1*(a4*(1-a4**3)*a2**3*(1-a2)
     &            -a2*(1-a2**3)*a4**3*(1-a4))
           a22=x1**2*(a4**2*(1-a4**2)*a2**3*(1-a2)
     &               -a2**2*(1-a2**2)*a4**3*(1-a4))

           b1=(y3-y1*a3**4)*a2**3*(1-a2)
     &       -(y2-y1*a2**4)*a3**3*(1-a3)
           b2=(y4-y1*a4**4)*a2**3*(1-a2)
     &       -(y2-y1*a2**4)*a4**3*(1-a4)

           c1=(w3-w1*a3**4)*a2**3*(1-a2)
     &       -(w2-w1*a2**4)*a3**3*(1-a3)
           c2=(w4-w1*a4**4)*a2**3*(1-a2)
     &       -(w2-w1*a2**4)*a4**3*(1-a4)

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dvxdy=(b1*a22-b2*a12)/det
           dvydy=(c1*a22-c2*a12)/det
       else
	write(6,*)'!!',det,a11,a22,a21,a12
           dvxdy=(vx(i,j+1)-vx(i,j-1))/dy(j)*0.5
           dvydy=(vy(i,j+1)-vy(i,j-1))/dy(j)*0.5
       endif

           x1= dx(i)/2+dx(i+1)/2
           x2= dx(i)/2+dx(i+1)+dx(i+2)/2
           x3=-dx(i)/2-dx(i-1)/2
           x4=-dx(i)/2-dx(i-1)-dx(i-2)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=vx(i+1,j)-vx(i,j)
           y2=vx(i+2,j)-vx(i,j)
           y3=vx(i-1,j)-vx(i,j)
           y4=vx(i-2,j)-vx(i,j)

           w1=vy(i+1,j)-vy(i,j)
           w2=vy(i+2,j)-vy(i,j)
           w3=vy(i-1,j)-vy(i,j)
           w4=vy(i-2,j)-vy(i,j)

           a11=x1*(a3*(1-a3**3)*a2**3*(1-a2)
     &            -a2*(1-a2**3)*a3**3*(1-a3))
           a12=x1**2*(a3**2*(1-a3**2)*a2**3*(1-a2)
     &               -a2**2*(1-a2**2)*a3**3*(1-a3))
           a21=x1*(a4*(1-a4**3)*a2**3*(1-a2)
     &            -a2*(1-a2**3)*a4**3*(1-a4))
           a22=x1**2*(a4**2*(1-a4**2)*a2**3*(1-a2)
     &               -a2**2*(1-a2**2)*a4**3*(1-a4))

           b1=(y3-y1*a3**4)*a2**3*(1-a2)
     &       -(y2-y1*a2**4)*a3**3*(1-a3)
           b2=(y4-y1*a4**4)*a2**3*(1-a2)
     &       -(y2-y1*a2**4)*a4**3*(1-a4)

           c1=(w3-w1*a3**4)*a2**3*(1-a2)
     &       -(w2-w1*a2**4)*a3**3*(1-a3)
           c2=(w4-w1*a4**4)*a2**3*(1-a2)
     &       -(w2-w1*a2**4)*a4**3*(1-a4)
          
       det=a11*a22-a21*a12

       if (abs(det).gt.1.d-20) then
           dvxdx=(b1*a22-b2*a12)/det
           dvydx=(c1*a22-c2*a12)/det
       else
        write(6,*)'xx',det,a11,a22,a12,a21
           dvxdx=(vx(i+1,j)-vx(i-1,j))/dx(i)*0.5
           dvydx=(vy(i+1,j)-vy(i-1,j))/dx(i)*0.5
       endif

         visxx(i,j)=visc(i,j)*(2.d0*dvxdx - (dvxdx + dvydy)*twoth)
         visxy(i,j)=visc(i,j)*(dvxdy + dvydx)
	 visyy(i,j)=visc(i,j)*(2.d0*dvydy - (dvxdx + dvydy)*twoth)


         enddo
         enddo

          return
          end


