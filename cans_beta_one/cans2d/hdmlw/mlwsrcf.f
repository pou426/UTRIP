c======================================================================|
      subroutine mlwsrcf(du,dt,s,ux0,ux1,ix,uy0,uy1,jx)
c======================================================================|
c
c NAME  mlwsrcf
c
c PURPOSE
c    apply source term for second half of Modified Lax-Wendroff method
c
c INPUTS & OUTPUTS
c    du(ix,jx): [double] variation in this step
c
c OUTPUTS
c    None
c
c INPUTS
c    s(ix,jx) : [double] source term
c    ux0(ix) : [double] 0.5*dx(i)/dxm(i)
c    ux1(ix) : [double] 0.5*dx(i-1)/dxm(i)
c    dt: [double] delta time 
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama based on R. Matsumoto's code
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension du(ix,jx),s(ix,jx)
      dimension ux0(ix),ux1(ix)
      dimension uy0(jx),uy1(jx)
c----------------------------------------------------------------------|      
      do j=2,jx-1
      do i=2,ix-1
         sh   = uy0(j)*(ux0(i)*s(i-1,j-1)+ux1(i)*s(i,j-1))
     &        + uy1(j)*(ux0(i)*s(i-1,j)  +ux1(i)*s(i,j)  )
c     summation of all terms
         du(i,j)= du(i,j)+0.5*dt*sh
      enddo
      enddo
c      
      return
      end
