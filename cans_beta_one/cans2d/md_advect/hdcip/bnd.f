c======================================================================|
      subroutine bnd(margin,ro,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix,jx)
      call bdsppy(0,margin,ro,ix,jx)
      call bdsppx(1,margin,ro,ix,jx)
      call bdsppy(1,margin,ro,ix,jx)

      return
      end
