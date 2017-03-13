c     -----------------------------------------------------------------
      subroutine resid(res,u,rhs,nx)
      integer nx
      double precision res(-1:nx,-1:nx),rhs(-1:nx,-1:nx),u(-1:nx,-1:nx)
      integer i,j
      double precision h,h2i
      h=1.d0/nx
      h2i=1.d0/(h*h)
c     --- boundary condition: u
c     --- begin: periodic boundary condition
      do i=-1,nx
         u(i,-1)=u(i,nx-1)
      enddo
      do j=-1,nx
         u(-1,j)=u(nx-1,j)
      enddo
c     --- end: periodic boundary condition
c
      do j=0,nx-1
         do i=0,nx-1
            res(i,j)=-h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-
     $           4.d0*u(i,j))+rhs(i,j)
         enddo
      enddo
c     --- boundary condition: res
c     --- begin: periodic boundary condition
      do i=0,nx
         res(i,nx)=res(i,0)
      enddo
      do j=0,nx
         res(nx,j)=res(0,j)
      enddo
c     --- end: periodic boundary condition
      return
      end
