c***********************************************************************
      subroutine bndcmat(margin,cmat,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7)
c======================================================================@

      if (ipe.eq.0) then
      call bmsppx(0,margin,cmat,ix,jx,kx)
      endif
      if (ipe.eq.ipex-1) then
      call bmsppx(1,margin,cmat,ix,jx,kx)
      endif

      if (jpe.eq.0) then
      call bmsppy(0,margin,cmat,ix,jx,kx)
      endif
      if (jpe.eq.jpex-1) then
      call bmsppy(1,margin,cmat,ix,jx,kx)
      endif

      if (kpe.eq.0) then
      call bmsppz(0,margin,cmat,ix,jx,kx)
      endif
      if (kpe.eq.kpex-1) then
      call bmsppz(1,margin,cmat,ix,jx,kx)
      endif

      return
      end
