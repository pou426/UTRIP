c======================================================================|
      subroutine bnd(margin,pr,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension pr(ix)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,pr,ix)

      call bdsppx(1,margin,pr,ix)

      return
      end
