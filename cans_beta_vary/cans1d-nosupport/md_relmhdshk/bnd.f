c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,by,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,pr,ix)
      call bdfrex(0,margin,ro,ix)
      call bdfrex(0,margin,vx,ix)
      call bdfrex(0,margin,vy,ix)
      call bdfrex(0,margin,by,ix)

      call bdfrex(1,margin,pr,ix)
      call bdfrex(1,margin,ro,ix)
      call bdfrex(1,margin,vx,ix)
      call bdfrex(1,margin,vy,ix)
      call bdfrex(1,margin,by,ix)

      return
      end
