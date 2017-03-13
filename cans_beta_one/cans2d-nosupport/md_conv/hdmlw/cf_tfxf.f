c======================================================================|
      subroutine cf_tfxf(tdfx,tdfy,te,dxm,ix,dym,jx,rkapm,teb1,tebx)
c======================================================================|
c
c NAME  cf_tfxf
c
c PURPOSE
c    calculate conduction flux
c    temperature is fixed at y boundaries
c    4th order inside and at y boundaries    
c
c    for the second step of LW scheme
c
c OUTPUTS
c    tdf?(ix,jx): [double] conduciton flux
c
c INPUTS
c    te(ix,jx): [double] temperature
c    dxm(ix),dym(jx) : [double] grid spacing (on i+1/2 grids)
c    ix,jx: [integer] dimension size
c    teb1, tebx ; [double] temperature at y boudaries
c
c HISTORY
c    written 2005/6/29 H.Isobe 
c
c----------------------------------------------------------------------|
      implicit real*8 (a-h,o-z)
      dimension dxm(ix),dym(jx)
      dimension te(ix,jx),rkapm(ix,jx)
      dimension tdfx(ix,jx),tdfy(ix,jx)
c----------------------------------------------------------------------|

         do j=3,jx-2
         do i=3,ix-2
           x1= dym(j)/2+dym(j+1)/2
           x2= dym(j)/2+dym(j+1)+dym(j+2)/2
           x3=-dym(j)/2-dym(j-1)/2
           x4=-dym(j)/2-dym(j-1)-dym(j-2)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=te(i,j+1)-te(i,j)
           y2=te(i,j+2)-te(i,j)
           y3=te(i,j-1)-te(i,j)
           y4=te(i,j-2)-te(i,j)

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

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dtedy=(b1*a22-b2*a12)/det
       else
           dtedy=(te(i,j+1)-te(i,j-1))/dym(j)*0.5
       endif

           x1= dxm(i)/2+dxm(i+1)/2
           x2= dxm(i)/2+dxm(i+1)+dxm(i+2)/2
           x3=-dxm(i)/2-dxm(i-1)/2
           x4=-dxm(i)/2-dxm(i-1)-dxm(i-2)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=te(i+1,j)-te(i,j)
           y2=te(i+2,j)-te(i,j)
           y3=te(i-1,j)-te(i,j)
           y4=te(i-2,j)-te(i,j)

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

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dtedx=(b1*a22-b2*a12)/det
       else
          dtedx=(te(i+1,j)-te(i-1,j))/dxm(i)*0.5d0
       endif
              
       tdfx(i,j)=rkapm(i,j)*dtedx
       tdfy(i,j)=rkapm(i,j)*dtedy

         enddo
         enddo

c----------------------------------------------------------------------|
c y boundaries (margin=2)
c 

         do i=1,ix

c --- j=3 ---
           j=3
           x1= dym(j)/2+dym(j+1)/2
           x2= dym(j)/2+dym(j+1)+dym(j+2)/2
           x3= dym(j)/2+dym(j+1)+dym(j+2)+dym(j+3)/2
           x4=-dym(j)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=te(i,j+1)-te(i,j)
           y2=te(i,j+2)-te(i,j)
           y3=te(i,j+3)-te(i,j)
           y4=teb1-te(i,j)

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

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dtedy3=(b1*a22-b2*a12)/det
       else
           dtedy3=(te(i,j+1)-te(i,j-1))/dym(j)*0.5d0
       endif

c --- j=4 ---

           j=4
           x1= dym(j)/2+dym(j+1)/2
           x2= dym(j)/2+dym(j+1)+dym(j+2)/2
           x3=-dym(j)/2-dym(j-1)/2
           x4=-dym(j)/2-dym(j-1)

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=te(i,j+1)-te(i,j)
           y2=te(i,j+2)-te(i,j)
           y3=te(i,j-1)-te(i,j)
           y4=teb1-te(i,j)

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

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dtedy4=(b1*a22-b2*a12)/det
       else
           dtedy4=(te(i,j+1)-te(i,j-1))/dym(j)*0.5d0
       endif

            tdfy(i,3)=rkapm(i,3)*dtedy3
            tdfy(i,4)=rkapm(i,4)*dtedy4


c --- j=jx-3 ---
           j=jx-3

           x1= dym(j)/2
           x2=-dym(j)/2-dym(j-1)/2
           x3=-dym(j)/2-dym(j-1)-dym(j-2)/2
           x4=-dym(j)/2-dym(j-1)-dym(j-2)-dym(j-3)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=tebx-te(i,j)
           y2=te(i,j-1)-te(i,j)
           y3=te(i,j-2)-te(i,j)
           y4=te(i,j-3)-te(i,j)

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

       det=a11*a22-a21*a12
       if (abs(det).gt.1.d-20) then
           dtedyx3=(b1*a22-b2*a12)/det
       else
           dtedyx3=(te(i,j+1)-te(i,j-1))/dym(j)*0.5d0
       endif

c --- j=jx-4 ---

           j=jx-4
           x1= dym(j)/2+dym(j+1)/2
           x2= dym(j)/2+dym(j+1)
           x3=-dym(j)/2-dym(j-1)/2
           x4=-dym(j)/2-dym(j-1)-dym(j-2)/2

           a2=x2/x1
           a3=x3/x1
           a4=x4/x1

           y1=te(i,j+1)-te(i,j)
           y2=tebx-te(i,j)
           y3=te(i,j-1)-te(i,j)
           y4=te(i,j-2)-te(i,j)

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

           det=a11*a22-a21*a12
           if (abs(det).gt.1.d-20) then
               dtedyx4=(b1*a22-b2*a12)/det
           else
               dtedyx4=(te(i,j+1)-te(i,j-1))/dym(j)*0.5d0
           endif

           tdfy(i,jx-3)=rkapm(i,jx-3)*dtedyx3
           tdfy(i,jx-4)=rkapm(i,jx-4)*dtedyx4


c outside the boundaries
            tdfy(i,1)=0.d0
            tdfy(i,2)=0.d0
            tdfx(i,1)=0.d0
            tdfx(i,2)=0.d0
            tdfy(i,jx-2)=0.d0
            tdfy(i,jx-1)=0.d0
            tdfy(i,jx)  =0.d0
            tdfx(i,jx-2)=0.d0
            tdfx(i,jx-1)=0.d0
            tdfx(i,jx)  =0.d0
         enddo


          return
          end
