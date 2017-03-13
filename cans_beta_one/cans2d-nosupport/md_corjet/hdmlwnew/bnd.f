c======================================================================|
      subroutine bnd(margin,ro,pr,vx,vy,bx,by,az,ix,jx,roi,pri)
c======================================================================|
c     apply boundary condition 
c----------------------------------------------------------------------|      
      implicit double precision (a-h,o-z)
      dimension ro(ix,jx),pr(ix,jx),vx(ix,jx),vy(ix,jx)
      dimension bx(ix,jx),by(ix,jx),az(ix,jx)
      dimension roi(ix,jx),pri(ix,jx)
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
      call bdsppx(1,margin,bx,ix,jx)
      call bdspnx(1,margin,by,ix,jx)
      call bdsppx(1,margin,az,ix,jx)

      call bdiniy(0,margin,ro,roi,ix,jx)
      call bdiniy(0,margin,pr,pri,ix,jx)
      call bdsppy(0,margin,vx,ix,jx)
      call bdspny(0,margin,vy,ix,jx)
      call bdsppy(0,margin,bx,ix,jx)
      call bdspny(0,margin,by,ix,jx)
      call bdsppy(0,margin,az,ix,jx)

      call bdiniy(1,margin,ro,roi,ix,jx)
      call bdiniy(1,margin,pr,pri,ix,jx)
      call bdsppy(1,margin,vx,ix,jx)
      call bdspny(1,margin,vy,ix,jx)
      call bdsppy(1,margin,bx,ix,jx)
      call bdspny(1,margin,by,ix,jx)
      call bdsppy(1,margin,az,ix,jx)

      return
      end
