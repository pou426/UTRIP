      subroutine ggtoggm(gx,gxm,gy,gym,ix,jx)

      implicit double precision (a-h,o-z)
      dimension gx(ix,jx),gy(ix,jx)
      dimension gxm(ix,jx),gym(ix,jx)

      do j=1,jx
      do i=1,ix-1
         gxm(i,j)=0.5*(gx(i,j)+gx(i+1,j))
      enddo
      enddo
      do j=1,jx-1
      do i=1,ix
         gym(i,j)=0.5*(gy(i,j)+gy(i,j+1))
      enddo
      enddo


      return
      end
