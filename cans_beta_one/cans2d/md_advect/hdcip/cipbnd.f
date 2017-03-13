c======================================================================|
      subroutine cipbnd(margin,ro,rodx,rody,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),rodx(ix,jx),rody(ix,jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,ro,ix,jx)
      call bdspnx(0,margin,rodx,ix,jx)
      call bdsppx(0,margin,rody,ix,jx)

      call bdsppx(1,margin,ro,ix,jx)
      call bdspnx(1,margin,rodx,ix,jx)
      call bdsppx(1,margin,rody,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,rodx,ix,jx)
      call bdspny(0,margin,rody,ix,jx)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,rodx,ix,jx)
      call bdspny(1,margin,rody,ix,jx)

      return
      end

