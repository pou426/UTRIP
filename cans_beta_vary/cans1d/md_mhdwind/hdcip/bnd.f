c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,by,gm,vstar,ix)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix),pr(ix),vx(ix),vy(ix),by(ix)
c----------------------------------------------------------------------|      
      ro0=1.d0
      pr0=1.d0/gm
      call bdcnsx(0,margin,ro,ro0,ix)
      call bdcnsx(0,margin,pr,pr0,ix)
      call bdfrex(0,margin,vx,ix)
      call bdcnsx(0,margin,vy,vstar,ix)
      call bdfrex(0,margin,by,ix)

      call bdfrex(1,margin,ro,ix)
      call bdfrex(1,margin,pr,ix)
      call bdfrex(1,margin,vx,ix)
      call bdfrex(1,margin,vy,ix)
      call bdfrex(1,margin,by,ix)

      return
      end
