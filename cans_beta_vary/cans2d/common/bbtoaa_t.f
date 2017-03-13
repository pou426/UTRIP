c======================================================================|
      subroutine bbtoaa_t(ay,bz,bx,dzm,dxm,margin,ix,jx)
c======================================================================|
c
c NAME  bbtoaa_t
c
c PURPOSE
c    calculate magnetic Vector potential
c
c OUTPUTS
c    ay(ix,jx): [double] magnetic Vector potential
c
c INPUTS
c    bz(ix,jx),bx(ix,jx): [double] magnetic field
c    dxm(ix),dzm(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dzm(jx)
      dimension bx(ix,jx),bz(ix,jx)
      dimension ay(ix,jx)
c----------------------------------------------------------------------|

      ay(margin+1,margin+1)=0.d0

      i=margin+1
      do j=margin+2,jx
        ay(i,j) = ay(i,j-1)
     &          - 0.5*(bx(i,j-1)+bx(i,j  ))*dzm(j-1)
      enddo
      do j=margin,1,-1
        ay(i,j) = ay(i,j+1)
     &          + 0.5*(bx(i,j  )+bx(i,j+1))*dzm(j-1)
      enddo

      do j=1,jx
      do i=margin+2,ix
         ay(i,j) = ay(i-1,j)
     &          + 0.5*(bz(i-1,j)+bz(i  ,j))*dxm(i-1)
      enddo
      do i=margin,1,-1
         ay(i,j) = ay(i+1,j)
     &          - 0.5*(bz(i  ,j)+bz(i+1,j))*dxm(i  )
      enddo
      enddo


          return
          end
