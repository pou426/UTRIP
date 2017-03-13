c***********************************************************************
      subroutine bndcmat(margin,cmat,ix,jx)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5)
c======================================================================@

      call bmsppx(0,margin,cmat,ix,jx)
      call bmsppx(1,margin,cmat,ix,jx)

      call bmsppy(0,margin,cmat,ix,jx)
      call bmsppy(1,margin,cmat,ix,jx)

      return
      end
