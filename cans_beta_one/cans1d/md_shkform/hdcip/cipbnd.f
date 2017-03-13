c======================================================================|
      subroutine cipbnd(margin,te,vxm,rodx,tedx,vxdxm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension te(ix),vxm(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix)
c----------------------------------------------------------------------|      
      call bdperx(margin,margin,rodx,ix)
      call bdperx(margin,margin,te,ix)
      call bdperx(margin,margin,tedx,ix)
      call bdperx(margin-1,margin,vxm,ix)
      call bdperx(margin-1,margin,vxdxm,ix)

      return
      end
