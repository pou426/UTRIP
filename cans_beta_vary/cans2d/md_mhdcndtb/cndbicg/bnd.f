c======================================================================|
      subroutine bnd(margin,pr,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pr(ix,jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,pr,ix,jx)

      call bdsppx(1,margin,pr,ix,jx)

      call bdsppy(0,margin,pr,ix,jx)

      call bdsppy(1,margin,pr,ix,jx)

      return
      end
