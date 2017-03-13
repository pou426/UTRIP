c***** self-gravity
      subroutine selfgper(ro,uod,gx,gy,x,y,ng,nx)
      implicit double precision (a-h,o-z)
      double precision x(-1:nx+1),y(-1:nx+1)
      double precision ro(-1:nx+1,-1:nx+1),
     1              gx(-1:nx+1,-1:nx+1),gy(-1:nx+1,-1:nx+1)
      integer ng,nx,memlen,ncycle,i,j
      double precision 
     &     u(-1:nx,-1:nx),ro7(-1:nx,-1:nx),
     &     uod(-1:nx,-1:nx)
      double precision h,h2,h2i
c     ---- definition of parameter ----
      h=1.0d0/nx
      h2=h*h
      h2i=1.d0/h2
      ncycle=2   !ncycle=2 -> w-cycle
      memlen=13*(4**ng-1)/3+28*2**ng+20*ng-55      
c     --- calculate mask
c     --- begin: periodic boundary condition
      uave=0.0d0
      do j=0,nx-1
         do i=0,nx-1
            uave=uave+ro(i,j)*h2
         enddo
      enddo
      do j=0,nx
         do i=0,nx
            ro7(i,j)=ro(i,j)-uave
         enddo
      enddo
c     --- boundary condition: ro7
      do i=0,nx
         ro7(i,-1)=ro7(i,nx-1)
         ro7(i,nx)=ro7(i,0)
      enddo
      do j=0,nx
         ro7(-1,j)=ro7(nx-1,j)
         ro7(nx,j)=ro7(0,j)
      enddo
         ro7(-1,-1)=ro7(nx-1,nx-1)
c     --- end: periodic boundary condition
c     ---- residual ----
      call resid(u,uod,ro7,nx)
c     ---- solve the pde ---
      call mglin(u,nx,ncycle,ng,memlen)
c     ---- update the potential & store the potential variables
      do j=0,nx
         do i=0,nx
            u(i,j)=u(i,j)+uod(i,j)
            uod(i,j)=u(i,j)
         enddo
      enddo
c     --- boundary condition: u
c     --- begin: periodic boundary condition
      do i=0,nx
         u(i,-1)=u(i,nx-1)
         u(i,nx)=u(i,0)
      enddo
      do j=0,nx
         u(-1,j)=u(nx-1,j)
         u(nx,j)=u(0,j)
      enddo
         u(-1,-1)=u(nx-1,nx-1)
c     --- end: periodic boundary condition
c     ---- calculation of gradient ----
      do j=0,nx-1
         do i=0,nx-1
            gx(i,j)=-(x(nx)-x(0))*(u(i+1,j)-u(i-1,j))/h
         end do
      end do
      do j=0,nx-1
         do i=0,nx-1
            gy(i,j)=-(y(nx)-y(0))*(u(i,j+1)-u(i,j-1))/h
         end do
      end do
c     --- boundary condition: g
c     --- begin: periodic boundary condition
      do i=0,nx-1
         gx(i,nx)=gx(i,0)
         gy(i,nx)=gy(i,0)
      enddo
      do j=0,nx
         gx(nx,j)=gx(0,j)
         gy(nx,j)=gy(0,j)
      enddo
c     --- end: periodic boundary condition
      return
      end
