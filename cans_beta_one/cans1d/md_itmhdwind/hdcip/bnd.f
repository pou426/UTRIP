c======================================================================|
      subroutine bnd(margin,ro,vx,vy,by,vstar,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),vx(ix),vy(ix),by(ix)
c----------------------------------------------------------------------|      
      ro0=1.d0
      call bdcnsx(0,margin,ro,ro0,ix)
      call bdfrex(0,margin,vx,ix)
      call bdcnsx(0,margin,vy,vstar,ix)
      call bdfrex(0,margin,by,ix)

      call bdfrex(1,margin,ro,ix)
      call bdfrex(1,margin,vx,ix)
      call bdfrex(1,margin,vy,ix)
      call bdfrex(1,margin,by,ix)

      return
      end
