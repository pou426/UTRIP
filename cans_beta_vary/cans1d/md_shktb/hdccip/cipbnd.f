c======================================================================|
      subroutine cipbnd(margin,te,vxm,rodx,tedx,vxdxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),te(ix),vxm(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix)
c----------------------------------------------------------------------|      
      call bdspnx(0,margin,rodx,ix)
      call bdsppx(0,margin,te,ix)
      call bdspnx(0,margin,tedx,ix)
      call bdsmnx(0,margin-1,vxm,ix)
      call bdsmpx(0,margin-1,vxdxm,ix)

      call bdspnx(1,margin,rodx,ix)
      call bdsppx(1,margin,te,ix)
      call bdspnx(1,margin,tedx,ix)
      call bdsmnx(1,margin,vxm,ix)
      call bdsmpx(1,margin,vxdxm,ix)

      return
      end
