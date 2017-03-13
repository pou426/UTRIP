c======================================================================|
      subroutine bnd(margin,ro,vx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix)
      dimension vx(ix)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix)
      call bdspnx(0,margin,vx,ix)

      call bdsppx(1,margin,ro,ix)
      call bdspnx(1,margin,vx,ix)


      return
      end
