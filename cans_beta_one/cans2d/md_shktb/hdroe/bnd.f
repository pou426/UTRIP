c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)
c----------------------------------------------------------------------|      
      call bdfrex(0,margin,ro,ix,jx)
      call bdfrex(0,margin,pr,ix,jx)
      call bdfrex(0,margin,vx,ix,jx)
      call bdfrex(0,margin,vy,ix,jx)

      call bdfrex(1,margin,ro,ix,jx)
      call bdfrex(1,margin,pr,ix,jx)
      call bdfrex(1,margin,vx,ix,jx)
      call bdfrex(1,margin,vy,ix,jx)

      call bdfrey(0,margin,ro,ix,jx)
      call bdfrey(0,margin,pr,ix,jx)
      call bdfrey(0,margin,vx,ix,jx)
      call bdfrey(0,margin,vy,ix,jx)

      call bdfrey(1,margin,ro,ix,jx)
      call bdfrey(1,margin,pr,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdfrey(1,margin,vy,ix,jx)

      return
      end
