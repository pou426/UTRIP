c***********************************************************************
      subroutine bndcmat(margin,cmat,ix,jx)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5)
c======================================================================@

      call bmperx(margin,margin,cmat,ix,jx)
      call bmpery(margin,margin,cmat,ix,jx)

      return
      end
