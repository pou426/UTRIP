c***********************************************************************
      subroutine bndcmat(margin,cmat,ix,jx,kx)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7)
c======================================================================@

      call bmsppx(0,margin,cmat,ix,jx,kx)
      call bmsppx(1,margin,cmat,ix,jx,kx)

      call bmsppy(0,margin,cmat,ix,jx,kx)
      call bmsppy(1,margin,cmat,ix,jx,kx)

      call bmsppz(0,margin,cmat,ix,jx,kx)
      call bmsppz(1,margin,cmat,ix,jx,kx)

      return
      end
