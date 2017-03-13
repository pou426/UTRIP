c======================================================================|
      subroutine bnd(margin,pr,ix,jx,kx
     &           ,ipe,jpe,kpe,ipex,jpex,kpex)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pr(ix,jx,kx)
c----------------------------------------------------------------------|      
      if (ipe.eq.0) then
      call bdsppx(0,margin,pr,ix,jx,kx)
      endif
      if (ipe.eq.ipex-1) then
      call bdsppx(1,margin,pr,ix,jx,kx)
      endif

      if (jpe.eq.0) then
      call bdsppy(0,margin,pr,ix,jx,kx)
      endif
      if (jpe.eq.jpex-1) then
      call bdsppy(1,margin,pr,ix,jx,kx)
      endif

      if (kpe.eq.0) then
      call bdsppz(0,margin,pr,ix,jx,kx)
      endif
      if (kpe.eq.kpex-1) then
      call bdsppz(1,margin,pr,ix,jx,kx)
      endif

      return
      end
