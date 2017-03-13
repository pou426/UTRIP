c======================================================================|
      subroutine gptogg(gx,gy,gxm,gym,gp,dx,dxm,dy,dym,ix,jx)
c======================================================================|
c
c NAME  gptogg
c
c PURPOSE
c    calculate gravity from potential
c
c OUTPUTS
c    gx(ix,jx),gy,gxm,gym: [double] gravitational acceleration
c
c INPUTS
c    gp(ix,jx): [double] gravitational potential
c    dx(ix),dy(jx),dxm(ix),dym(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2005-2-21 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dx(ix),dy(jx),dxm(ix),dym(jx)
      dimension gx(ix,jx),gy(ix,jx),gxm(ix,jx),gym(ix,jx)
      dimension gp(ix,jx)
c----------------------------------------------------------------------|

      do i=2,ix-1
      do j=2,jx-1
        gx(i,j)=-(gp(i+1,j)-gp(i-1,j))/dx(i)
        gy(i,j)=-(gp(i,j+1)-gp(i,j-1))/dy(j)
      enddo
      enddo

      do i=1,ix-1
      do j=1,jx-1
        gxm(i,j)=-(gp(i+1,j)-gp(i,j))/dxm(i)
        gym(i,j)=-(gp(i,j+1)-gp(i,j))/dym(j)
      enddo
      enddo


      return
      end
