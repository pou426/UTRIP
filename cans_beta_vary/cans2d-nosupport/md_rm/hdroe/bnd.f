c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx)
      call bdfrex(0,margin,pr,ix,jx)
      call bdfrex(0,margin,vx,ix,jx)
      call bdfrex(0,margin,vy,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)

      call bdpery(margin,margin,ro,ix,jx)
      call bdpery(margin,margin,pr,ix,jx)
      call bdpery(margin,margin,vx,ix,jx)
      call bdpery(margin,margin,vy,ix,jx)

      return
      end

