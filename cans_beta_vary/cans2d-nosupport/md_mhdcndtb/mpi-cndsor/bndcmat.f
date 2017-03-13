c***********************************************************************
      subroutine bndcmat(margin,cmat,ix,jx
     &           ,ipe,jpe,ipex,jpex)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5)
c======================================================================@

      if (ipe.eq.0) then
      call bmsppx(0,margin,cmat,ix,jx)
      endif
      if (ipe.eq.ipex-1) then
      call bmsppx(1,margin,cmat,ix,jx)
      endif

      if (jpe.eq.0) then
      call bmsppy(0,margin,cmat,ix,jx)
      endif
      if (jpe.eq.jpex-1) then
      call bmsppy(1,margin,cmat,ix,jx)
      endif


      return
      end
