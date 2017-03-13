c======================================================================|
      subroutine cipbnd(margin,rodx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension rodx(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,rodx,ix)

      return
      end

