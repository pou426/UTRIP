c======================================================================|
      subroutine bnd(margin,pr,ix,jx,kx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pr(ix,jx,kx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,pr,ix,jx,kx)
      call bdsppx(1,margin,pr,ix,jx,kx)

      call bdsppy(0,margin,pr,ix,jx,kx)
      call bdsppy(1,margin,pr,ix,jx,kx)

      call bdsppz(0,margin,pr,ix,jx,kx)
      call bdsppz(1,margin,pr,ix,jx,kx)

      return
      end
