c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,bx,by,az,ix,jx,dym)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx),az(ix,jx)
      dimension dym(jx)
c----------------------------------------------------------------------|      
      call bdsppx(0,margin,pr,ix,jx)
      call bdsppx(0,margin,ro,ix,jx)
      call bdspnx(0,margin,vx,ix,jx)
      call bdsppx(0,margin,vy,ix,jx)
      call bdspnx(0,margin,bx,ix,jx)
      call bdsppx(0,margin,by,ix,jx)
      call bdsppx(0,margin,az,ix,jx)

      call bdsppx(1,margin,pr,ix,jx)
      call bdsppx(1,margin,ro,ix,jx)
      call bdspnx(1,margin,vx,ix,jx)
      call bdsppx(1,margin,vy,ix,jx)
      call bdspnx(1,margin,bx,ix,jx)
      call bdsppx(1,margin,by,ix,jx)
      call bdsppx(1,margin,az,ix,jx)

      call bdfrdy(0,margin,ro,dym,ix,jx)
      call bdfrdy(0,margin,pr,dym,ix,jx)
      call bdfrey(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)
      call bdfrey(0,margin,bx,ix,jx)
      call bdfrey(0,margin,by,ix,jx)
      call bdfrey(0,margin,az,ix,jx)

      call bdfrdy(1,margin,ro,dym,ix,jx)
      call bdfrdy(1,margin,pr,dym,ix,jx)
      call bdfrey(1,margin,vx,ix,jx)
      call bdspny(1,margin,vy,ix,jx)
      call bdfrey(1,margin,bx,ix,jx)
      call bdfrey(1,margin,by,ix,jx)
      call bdfrey(1,margin,az,ix,jx)

      return
      end
