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
      call bdperx(margin,margin,ro,ix,jx)
      call bdperx(margin,margin,pr,ix,jx)
      call bdperx(margin,margin,vx,ix,jx)
      call bdperx(margin,margin,vy,ix,jx)

      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)

      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,pr,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vy,ix,jx)

      return
      end
