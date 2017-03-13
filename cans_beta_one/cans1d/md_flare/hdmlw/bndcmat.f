c***********************************************************************
      subroutine bndcmat(margin,cmat,ix)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,3)
c======================================================================@

      call bmsppx(0,margin,cmat,ix)
      call bmsppx(1,margin,cmat,ix)

         return
         end
