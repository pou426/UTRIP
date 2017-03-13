      subroutine tvdminmod1(mdir,flux,fgrid,ix,jx)
      implicit double precision (a-h,o-z)

      dimension flux(ix,jx)
      dimension fgrid(ix,jx)
      dimension ftmp(ix,jx)

c----------------------------------------------------------------------|
c     define limiter functions
c     1. minmod limiter

      flmt(a,b)=max(0.0d0,min(b*sign(1.0d0,a),abs(a)))*sign(1.0d0,a)
c----------------------------------------------------------------------|

        do j=2,jx-2
        do i=2,ix-2
          ftmp(i,j)=flux(i,j)
        enddo
        enddo

      if (mdir.eq.1) then
        do j=2,jx-2
        do i=2,ix-2
          flux(i,j)=ftmp(i,j)
     &         +0.5*flmt(fgrid(i,j)-ftmp(i-1,j),fgrid(i+1,j)-ftmp(i,j))
     &         -0.5*flmt(ftmp(i+1,j)-fgrid(i+1,j),ftmp(i,j)-fgrid(i,j))
        enddo
        enddo
      endif

      if (mdir.eq.2) then
        do j=2,jx-2
        do i=2,ix-2
          flux(i,j)=ftmp(i,j)
     &         +0.5*flmt(fgrid(i,j)-ftmp(i,j-1),fgrid(i,j+1)-ftmp(i,j))
     &         -0.5*flmt(ftmp(i,j+1)-fgrid(i,j+1),ftmp(i,j)-fgrid(i,j))
        enddo
        enddo
      endif

      return
      end
