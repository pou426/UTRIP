c***********************************************************************
      subroutine bndcmat(cmat,ix,jx,kx,mfdim,margar)
c***********************************************************************
      implicit double precision (a-h,o-z)
      dimension mfdim(3),margar(3)
      dimension cmat(ix,jx,kx,7)
c======================================================================@

      if (mfdim(1).eq.1) then
      call bmsppx(0,margar(1),cmat,ix,jx,kx)
      call bmsppx(1,margar(1),cmat,ix,jx,kx)
      endif

      if (mfdim(2).eq.1) then
      call bmsppy(0,margar(2),cmat,ix,jx,kx)
      call bmsppy(1,margar(2),cmat,ix,jx,kx)
      endif

      if (mfdim(3).eq.1) then
      call bmsppz(0,margar(3),cmat,ix,jx,kx)
      call bmsppz(1,margar(3),cmat,ix,jx,kx)
      endif

      return
      end
