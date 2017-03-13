c======================================================================|
      subroutine cipbnd(margin,ro,vxm,rodx,vxdxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),vxm(ix)
      dimension rodx(ix),vxdxm(ix)
c----------------------------------------------------------------------|      
      zero=0.0d0

      call bdsppx(0,margin,ro,ix)
      call bdspnx(0,margin,rodx,ix)
      call bdsmnx(0,margin-1,vxm,ix)
      call bdsmpx(0,margin-1,vxdxm,ix)

      call bdsppx(1,margin,ro,ix)
      call bdspnx(1,margin,rodx,ix)
      call bdsmnx(1,margin,vxm,ix)
      call bdsmpx(1,margin,vxdxm,ix)

c     call bdfrex(1,margin,ro,ix)
c     call bdcnsx(1,margin,rodx,zero,ix)
c     call bdcnsx(1,margin,vxm,zero,ix)
c     call bdfrex(1,margin,vxdxm,ix)

      return
      end
