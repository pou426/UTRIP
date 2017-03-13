c======================================================================|
      subroutine cipbnd(margin,te,vxm,rodx,tedx,vxdxm,gm,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),te(ix),vxm(ix)
      dimension rodx(ix),tedx(ix),vxdxm(ix)
c----------------------------------------------------------------------|      
      zero=0.0d0
      te0=1.d0
      call bdcnsx(0,margin,rodx,zero,ix)
      call bdcnsx(0,margin,te,te0,ix)
      call bdcnsx(0,margin,tedx,zero,ix)
c
c  Don't comment out !
c  If following lines are commented out, then, 
c  artificial oscillation appears.
c
c     call bdcnsx(0,margin-1,vxm,zero,ix)
c     call bdcnsx(0,margin-1,vxdxm,zero,ix)

      call bdcnsx(1,margin,rodx,zero,ix)
      call bdfrex(1,margin,te,ix)
      call bdcnsx(1,margin,tedx,zero,ix)
      call bdfrex(1,margin,vxm,ix)
      call bdcnsx(1,margin,vxdxm,zero,ix)

      return
      end
