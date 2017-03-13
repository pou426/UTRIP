c======================================================================|
      subroutine mlwreset(u2,uh,u,ix,jx)
c======================================================================|
c
c NAME  mlwhalf
c
c PURPOSE
c    first half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c
c OUTPUTS
c    un(ix,jx) : [double] half step results on mid-grid points
c
c INPUTS
c    u(ix,jx) : [double] basic variables    
c    f(ix,jx) : [double] flux in x-direction
c    dxi(ix), dxim(ix) : [double] 1/dx
c    dt: [double] delta time 
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension u(ix,jx), uh(ix,jx), u2(ix,jx)

      do j=1,jx
      do i=1,ix
        u2(i,j)=u(i,j)
      enddo
      enddo

      do j=1,jx-1
      do i=1,ix-1
         uh(i,j)= (u(i+1,j)+u(i,j)+u(i+1,j+1)+u(i,j+1))/4.d0
      enddo
      enddo

      return
      end
