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

      call bdfrex(0,margin,ro,ix)
      call bdcnsx(0,margin,rodx,zero,ix)
      call bdfrex(0,margin-1,vxm,ix)
      call bdcnsx(0,margin-1,vxdxm,zero,ix)

      call bdfrex(1,margin,ro,ix)
      call bdcnsx(1,margin,rodx,zero,ix)
      call bdfrex(1,margin,vxm,ix)
      call bdcnsx(1,margin,vxdxm,zero,ix)

      return
      end
