c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,bx,by,az,ix,jx)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx)
      dimension pr(ix,jx)
      dimension vx(ix,jx)
      dimension vy(ix,jx)
      dimension bx(ix,jx)
      dimension by(ix,jx)
      dimension az(ix,jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,pr,ix,jx)
      call bdsppx(0,margin,ro,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdsppx(0,margin,vy,ix,jx)
      call bdsppx(0,margin,bx,ix,jx)
      call bdspnx(0,margin,by,ix,jx)
      call bdsppx(0,margin,az,ix,jx)

      call bdsppx(1,margin,pr,ix,jx)
      call bdsppx(1,margin,ro,ix,jx)
      call bdspnx(1,margin,vx,ix,jx)
      call bdsppx(1,margin,vy,ix,jx)
      call bdspnx(1,margin,bx,ix,jx)
      call bdsppx(1,margin,by,ix,jx)
      call bdsppx(1,margin,az,ix,jx)

      call bdsppy(0,margin,pr,ix,jx)
      call bdsppy(0,margin,ro,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)
      call bdspny(0,margin,bx,ix,jx)
      call bdsppy(0,margin,by,ix,jx)
      call bdsppy(0,margin,az,ix,jx)

      call bdsppy(1,margin,pr,ix,jx)
      call bdsppy(1,margin,ro,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vy,ix,jx)
      call bdspny(1,margin,bx,ix,jx)
      call bdsppy(1,margin,by,ix,jx)
      call bdsppy(1,margin,az,ix,jx)

      return
      end
