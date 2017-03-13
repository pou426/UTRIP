c======================================================================|
      subroutine cipbnd(margin,vxm,vy,by,rodx,vxdxm,vydx,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension vxm(ix),vy(ix),by(ix)
      dimension rodx(ix),vxdxm(ix),vydx(ix)
c----------------------------------------------------------------------|      
      call bdspnx(0,margin,rodx,ix)
      call bdsmnx(0,margin-1,vxm,ix)
      call bdsmpx(0,margin-1,vxdxm,ix)
      call bdsppx(0,margin,vy,ix)
      call bdspnx(0,margin,vydx,ix)
      call bdsppx(0,margin,by,ix)

      call bdspnx(1,margin,rodx,ix)
      call bdsmnx(1,margin,vxm,ix)
      call bdsmpx(1,margin,vxdxm,ix)
      call bdsppx(1,margin,vy,ix)
      call bdspnx(1,margin,vydx,ix)
      call bdsppx(1,margin,by,ix)

      return
      end
