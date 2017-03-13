c======================================================================|
      subroutine bnd(pr,ix,jx,kx,mfdim,margar)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension mfdim(3),margar(3)
      dimension pr(ix,jx,kx)
c----------------------------------------------------------------------|      
      if (mfdim(1).eq.1) then
      call bdsppx(0,margar(1),pr,ix,jx,kx)
      call bdsppx(1,margar(1),pr,ix,jx,kx)
      endif

      if (mfdim(2).eq.1) then
      call bdsppy(0,margar(2),pr,ix,jx,kx)
      call bdsppy(1,margar(2),pr,ix,jx,kx)
      endif

      if (mfdim(3).eq.1) then
      call bdsppz(0,margar(3),pr,ix,jx,kx)
      call bdsppz(1,margar(3),pr,ix,jx,kx)
      endif

      return
      end
