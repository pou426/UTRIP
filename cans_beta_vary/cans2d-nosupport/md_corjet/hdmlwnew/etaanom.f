c***********************************************************************
      subroutine etaanom(et,etm,bx,by,ro,dx,ix,dy,jx)
c***********************************************************************
c
c***********************************************************************
      implicit double precision (a-h,o-z)

      dimension et(ix,jx),etm(ix,jx)
      dimension cz(ix,jx)
      dimension ro(ix,jx),bx(ix,jx),by(ix,jx)
      dimension dx(ix),dy(jx)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      call bbtocz(cz,bx,by,dx,dy,ix,jx)

      etst=0.d0
      eta0=0.01d0
      vdcri=1000.d0
      etmax=1.d0
      rolim=1.d-4

      do j=1,jx
      do i=1,ix
        et(i,j)=0.d0
        etm(i,j)=0.d0
      enddo
      enddo

      do j=2,jx-1
      do i=2,ix-1
        if (ro(i,j).le.rolim) then
          vd = abs(cz(i,j)/ro(i,j))
          if(vd.gt.vdcri) et(i,j)=etst+eta0*(vd/vdcri-1.)**2
          if(et(i,j).ge.etmax) et(i,j)=etmax
        endif
      enddo
      enddo

      do j=2,jx-2
      do i=2,ix-2
        etm(i,j)=(et(i,j)+et(i,j+1)+et(i+1,j)+et(i+1,j+1))/4
      enddo
      enddo


      return
      end
