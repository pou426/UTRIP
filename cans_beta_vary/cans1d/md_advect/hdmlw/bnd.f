c======================================================================|
      subroutine bnd(margin,ro,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,ro,ix)

      return
      end
