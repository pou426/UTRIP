c======================================================================|
      subroutine cipbnd(margin,vxm,rodx,vxdxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension vxm(ix)
      dimension rodx(ix),vxdxm(ix)
c----------------------------------------------------------------------|      
      call bdspnx(0,margin,rodx,ix)
      call bdsmnx(0,margin-1,vxm,ix)
      call bdsmpx(0,margin-1,vxdxm,ix)

      call bdspnx(1,margin,rodx,ix)
      call bdsmnx(1,margin,vxm,ix)
      call bdsmpx(1,margin,vxdxm,ix)


      return
      end
