c     ----------------------------------------------------------------
      subroutine relax(u,rhs,nx)
      integer nx
      double precision rhs(-1:nx,-1:nx),u(-1:nx,-1:nx)
      integer i,j,ipass ,isw
      double precision h,h2
      h=1.d0/nx
      h2=h*h
      do ipass=1,2
c     --- boundary condition: u
c     --- begin: periodic boundary condition
         isw=mod(ipass,2)
         do i=isw,nx,2
            u(i,-1)=u(i,nx-1)
         enddo
         do j=isw,nx,2
            u(-1,j)=u(nx-1,j)
         enddo
c     --- end: periodic boundary condition
c     ---- main ----
         do j=0,nx-1
            isw=mod(j+ipass,2)
            do i=isw,nx-1,2
               u(i,j)=0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
     $              -h2*rhs(i,j))
            enddo
         enddo
c     --- boundary condition: u
c     --- begin: periodic boundary condition
         isw=mod(ipass,2)
         do i=isw,nx,2
            u(i,nx)=u(i,0)
         enddo
         do j=isw,nx,2
            u(nx,j)=u(0,j)
         enddo         
      enddo         
c     --- end: periodic boundary condition
      return
      end
