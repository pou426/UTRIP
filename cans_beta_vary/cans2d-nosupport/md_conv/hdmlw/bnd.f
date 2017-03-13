c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,ix,jx,gasr,teb1,tebx)
c======================================================================|
c     apply boundary condition 
c     symmetric at the boundary so that vy=0 and dvx/dy = 0
c     T is fixed at the boundary. 
c     Density is symmetric for outer region and solved by LW at 
c     y boundaries. 
c----------------------------------------------------------------------|      
      implicit real*8 (a-h,o-z)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)

c----------------------------------------------------------------------|      
      marginx=5

      do i=1,ix
         vx(i,margin+1)=vx(i,margin+2)
         vy(i,margin+1)=0.d0      
         pr(i,margin+1)=ro(i,margin+1)*teb1*gasr
         
         vx(i,jx-margin)=vx(i,jx-margin-1)
         vy(i,jx-margin)=0.d0
         pr(i,jx-margin)=ro(i,jx-margin)*tebx*gasr
   
      do j=1,margin
         vx(i,j)= vx(i,2*(margin+1)-j)
   	 vy(i,j)= 0.d0
         ro(i,j)= ro(i,2*(margin+1)-j)
         pr(i,j)= ro(i,j)*teb1*gasr

         vx(i,jx-margin+j)= vx(i,jx-margin-j)
	 vy(i,jx-margin+j)= 0.d0
         ro(i,jx-margin+j)= ro(i,jx-margin-j)
         pr(i,jx-margin+j)= ro(i,jx-margin+j)*tebx*gasr
      enddo

      enddo
            
      call bdperx(marginx,marginx,ro,ix,jx)
      call bdperx(marginx,marginx,pr,ix,jx)
      call bdperx(marginx,marginx,vx,ix,jx)
      call bdperx(marginx,marginx,vy,ix,jx)

      return
      end








